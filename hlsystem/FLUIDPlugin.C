#include <UT/UT_DSOVersion.h>
#include <UT/UT_Math.h>
#include <UT/UT_Interrupt.h>
#include <GU/GU_Detail.h>
#include <GA/GA_Primitive.h>
#include <GU/GU_PrimPoly.h>
#include <CH/CH_LocalVariable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_SpareData.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>

#include <GU/GU_PrimSphere.h>
#include <CH/CH_Manager.h>
#include <OP/OP_Director.h>
#include <OP/OP_AutoLockInputs.h>

#include <limits.h>
#include "FLUIDPlugin.h"

#include <HOM/HOM_ui.h>
#include <glm/gtx/string_cast.hpp>

#define FRAME_RANGE 100

void newSopOperator(OP_OperatorTable *table) {
    table->addOperator(
	    new OP_Operator("H2O-plugin",			// Internal name
			    "H2O",			// UI name
			     SOP_Fluid::myConstructor,	// How to build the SOP
			     SOP_Fluid::myTemplateList,	// My parameters
			     1,				// Min # of sources
			     1,				// Max # of sources
			     0,	// Local variables e.g., CH_LocalVariable SOP_Fluid::myVariables[] = {};
			     OP_FLAG_GENERATOR)		// Flag it as generator
	    );
}

static PRM_Name		constraintIteration("constraintIteration", "Constraint Iteration");
static PRM_Name		artificialPressure("artificialPressure", "Artificial Pressure");
static PRM_Name		viscosity("viscosity", "Viscosity");
static PRM_Name		vorticityConfinement("vorticityConfinement", "Vorticity Confinement");
//				     ^^^^^^^^    ^^^^^^^^^^^^^^^
//				     internal    descriptive version

// initial/default values for your parameters here
static PRM_Default constraintIterationDefault(2);
static PRM_Default artificialPressureDefault(0.0001);
static PRM_Default viscosityDefault(0.01);
static PRM_Default vorticityConfinementDefault(0.0000);

static PRM_Range iterationRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 30);
static PRM_Range tensileRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.01);
static PRM_Range viscosityRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.1);
static PRM_Range vorticityRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.001);

PRM_Template
SOP_Fluid::myTemplateList[] = {
	// default vals
	PRM_Template(PRM_INT,	PRM_Template::PRM_EXPORT_MIN, 1, &constraintIteration, &constraintIterationDefault, 0, &iterationRange),
	PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &artificialPressure, &artificialPressureDefault, 0, &tensileRange),
	PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &viscosity, &viscosityDefault, 0, &viscosityRange),
	PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &vorticityConfinement, &vorticityConfinementDefault, 0, &vorticityRange),
	PRM_Template()
};

// --------------------------end boilerplates-----------------------------------

OP_Node* SOP_Fluid::myConstructor(OP_Network *net, const char *name, OP_Operator *op) {
    return new SOP_Fluid(net, name, op);
}

SOP_Fluid::SOP_Fluid(OP_Network *net, const char *name, OP_Operator *op) : SOP_Node(net, name, op) {
	myFS = new FluidSystem();
	// SET SPH RAD - DUE to sensitivity of SPH sim, we REQUIRE 0.5 distance between points.
	myFS->SPH_RADIUS = 0.1;

	oldIteration = 2;
	oldKCorr = 0.0001;
	oldViscosity = 0.01;
	oldVorticity = 0.0003;

	init = true;
}

void SOP_Fluid::runSimulation(int frameNumber) {
	int moreFrame = frameNumber - totalPos.size();
	if (moreFrame > 0) {
		for (int i = 0; i <= moreFrame + FRAME_RANGE; ++i) {
			std::vector<glm::dvec3> temp;
			for (auto& f : myFS->fluidPs) {
				glm::dvec3 scaledPos = f->pos;
				scaledPos /= myFS->SPH_RADIUS;
				temp.push_back(scaledPos);
			}
			totalPos.push_back(temp);

			myFS->Run();
		}
	}
}

SOP_Fluid::~SOP_Fluid() {}

OP_ERROR SOP_Fluid::cookMySop(OP_Context &context) {
	OP_Node::flags().setTimeDep(true);// indicate that we have to cook every time current frame changes).
	fpreal now = context.getTime();
	fpreal currframe = OPgetDirector()->getChannelManager()->getSample(now);

	// update parameters
	int ite = CONSTRAINT_ITERATION(now);
	float kCorr = ARTIFICIAL_PRESSURE(now);
	float visc = VISCOSITY(now);
	float vorticity = VORTICITY_CONFINEMENT(now);


	//get iinputs

	//OP_AutoLockInputs inputs(this);
	//if (inputs.lockInput(0, context) >= UT_ERROR_ABORT)
	//	return error();
	if (lockInput(0, context) >= UT_ERROR_ABORT) { // check for 1 input (presumably geom)
		return error();
	}
	duplicateSource(0, context); // copy from input geometry to sop's own gdp
	// 1st connected input to node
	GU_Detail* fluid_gdp = new GU_Detail(inputGeo(0, context));

	//TODO: set magic number in ui? max particles
	if (int npts = fluid_gdp->getPointRange().getEntries() > 1000) {
		addWarning(SOP_MESSAGE, "Too many pts: " + npts);
		//HOM_Module& hou = HOM();
		//hou.ui().displayMessage("too many pts: " + npts);
		return error();
	}

	// inefficient to do this eveery time?
	std::vector<glm::dvec3> posn;
	GA_Offset ptoff;
	GA_FOR_ALL_PTOFF(fluid_gdp, ptoff) {
		UT_Vector3 pos = fluid_gdp->getPos3(ptoff);
		posn.push_back(glm::dvec3(pos[0], pos[2], pos[1]));

		// check points from volume dist? - somehow automate SPH_RAD?
		//double pDist = glm::length(posn.at(0) - posn.at(1));
	}
	if (init) { // do this just once in the beginning
		myFS->SPH_CreateExample(posn);
		init = false;
	}

	if (ite != oldIteration || kCorr != oldKCorr || visc != oldViscosity || vorticity != oldVorticity) {
		oldIteration = ite;
		oldKCorr = kCorr;
		oldViscosity = visc;
		oldVorticity = vorticity;
		myFS->setParameters(ite, visc, vorticity, kCorr);
		myFS->SPH_CreateExample(posn);
		totalPos.clear();
	}

	runSimulation(currframe);

	/*std::cout << "iteration " << ite << endl;
	cout << "stiffnes " << stif << endl;
	cout << "viscosity " << visc << endl;
	cout << "vorticity " << vorticity << "\n" << endl;*/

    UT_Vector4 pos;
    GU_PrimPoly		*poly;
    UT_Interrupt	*boss;

    if (error() < UT_ERROR_ABORT) {
		boss = UTgetInterrupt();
		gdp->clearAndDestroy();

		if (boss->opStart("Building Fluid") && 
			currframe < totalPos.size()) {	// currframe generation might not be able to catch up?
			for (auto& f : totalPos.at(currframe)) {
				glm::dvec3 scaledPos = f;
				//scaledPos /= SPH_RADIUS;

				UT_Vector3 pos;
				pos(0) = scaledPos.x;
				pos(1) = scaledPos.z;
				pos(2) = scaledPos.y;

				GU_PrimSphereParms sphere(gdp);
				double scale = myFS->SPH_RADIUS * 0.4;
				sphere.xform.scale(scale, scale, scale);
				sphere.xform.translate(pos);
				GU_PrimSphere::build(sphere, GEO_PRIMSPHERE);
			}
			// Highlight the star which we have just generated.  This routine
			// call clears any currently highlighted geometry, and then it
			// highlights every primitive for this SOP. 
			select(GU_SPrimitive);
		}
		// Tell the interrupt server that we've completed. Must do this
		// regardless of what opStart() returns.
		boss->opEnd();
    }

    return error();
}