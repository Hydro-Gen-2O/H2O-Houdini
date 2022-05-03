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

#define FRAME_RANGE 60 // a buffer for more frame than current frame

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
static PRM_Name		minCorner("minCorner", "Min Corner");
static PRM_Name		maxCorner("maxCorner", "Max Corner");
//static PRM_Name		maxPts("maxPts", "Maximum Points");
static PRM_Name		simulateButton("simulateButton", "Run Simulation");
//				     ^^^^^^^^    ^^^^^^^^^^^^^^^
//				     internal    descriptive version

// initial/default values for your parameters here
static PRM_Default constraintIterationDefault(2);
static PRM_Default artificialPressureDefault(0.0001);
static PRM_Default viscosityDefault(0.01);
static PRM_Default vorticityConfinementDefault(0.0000);
static PRM_Default minDefault[] = { PRM_Default(-10.0), PRM_Default(0.0), PRM_Default(-10.0) };
static PRM_Default maxDefault[] = { PRM_Default(10.0), PRM_Default(30.0), PRM_Default(10.0) };
//static PRM_Default maxPtsDefault(5000);

static PRM_Range iterationRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 30);
static PRM_Range tensileRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.01);
static PRM_Range viscosityRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.1);
static PRM_Range vorticityRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.001);
//static PRM_Range maxPtsRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 100000);

PRM_Template
SOP_Fluid::myTemplateList[] = {
	// default vals
	PRM_Template(PRM_INT,	PRM_Template::PRM_EXPORT_MIN, 1, &constraintIteration, &constraintIterationDefault, 0, &iterationRange),
	PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &artificialPressure, &artificialPressureDefault, 0, &tensileRange),
	PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &viscosity, &viscosityDefault, 0, &viscosityRange),
	PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &vorticityConfinement, &vorticityConfinementDefault, 0, &vorticityRange),
	PRM_Template(PRM_XYZ_J, 3, &minCorner, minDefault),
	PRM_Template(PRM_XYZ_J, 3, &maxCorner, maxDefault),
	//PRM_Template(PRM_INT,	PRM_Template::PRM_EXPORT_MIN, 1, &maxPts, &maxPtsDefault, 0, &maxPtsRange),
	PRM_Template(PRM_CALLBACK, 1, &simulateButton, 0, 0, 0, &simulate),
	PRM_Template()
};

int
SOP_Fluid::simulate(void* op, int index, fpreal t, const PRM_Template*)
{
	SOP_Fluid* fluid = (SOP_Fluid*)op;

	if (fluid->clearOrNot)
	{
		fluid->totalPos.clear();
		fluid->clearOrNot = false;
	}
	//fluid->totalPos.clear(); // clear cache
	fluid->runSimulation(fluid->currentFrame);
	//std::cout << fluid->totalPos.size() << std::endl;
	//fluid->cookMySop(*(fluid->myContext));
	//fluid->buildGeo();
	return 1;
}

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
	clearOrNot = false;

	init = true;
}

void SOP_Fluid::runSimulation(int frameNumber) {
	int moreFrame = frameNumber - totalPos.size();
	if (moreFrame >= 0) {
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
	currentFrame = currframe;
	myContext = &context;

	// update parameters
	int ite = CONSTRAINT_ITERATION(now);
	float kCorr = ARTIFICIAL_PRESSURE(now);
	float visc = VISCOSITY(now);
	float vorticity = VORTICITY_CONFINEMENT(now);
	float minx; //current mincorner
	float miny;
	float minz;
	minx = evalFloat("minCorner", 0, now);
	miny = evalFloat("minCorner", 1, now);
	minz = evalFloat("minCorner", 2, now);
	float maxx; // currentmaxcorner
	float maxy;
	float maxz;
	maxx = evalFloat("maxCorner", 0, now);
	maxy = evalFloat("maxCorner", 1, now);
	maxz = evalFloat("maxCorner", 2, now);
	glm::dvec3 currentMinCorner(minx, minz, miny);
	glm::dvec3 currentMaxCorner(maxx, maxz, maxy);


	//int maxPts = MAX_PTS(now);

	//get inputs
	OP_AutoLockInputs inputs(this);
	if (inputs.lockInput(0, context) >= UT_ERROR_ABORT) { return error(); }

	// only if input geo is different, re-get all the pts
	int input_changed;
	duplicateChangedSource(0, context, &input_changed);
	if (input_changed) {
		fluidPs.clear();
		// 1st connected input to node
		GU_Detail* fluid_gdp = new GU_Detail(inputGeo(0, context));
		/*if (int npts = fluid_gdp->getPointRange().getEntries() > maxPts) {
			addWarning(SOP_MESSAGE, "Too many pts: " + npts);
			return error();
		}*/
		
		GA_Offset ptoff;
		GA_FOR_ALL_PTOFF(fluid_gdp, ptoff) {
			UT_Vector3 pos = fluid_gdp->getPos3(ptoff);
			fluidPs.push_back(glm::dvec3(pos[0], pos[2], pos[1]));
			// check points from volume dist? - somehow automate SPH_RAD?
			//double pDist = glm::length(fluidPs.at(0) - fluidPs.at(1));
		}
		myFS->setParameters(ite, visc, vorticity, kCorr);
		myFS->SPH_CreateExample(fluidPs);
		clearOrNot = true;
	}
	
	if (init) { // do this just once in the beginning
		myFS->SPH_CreateExample(fluidPs);
		init = false;
	}

	if (ite != oldIteration || kCorr != oldKCorr || visc != oldViscosity || vorticity != oldVorticity ||
				currentMinCorner != myFS->SPH_VOLMIN || currentMaxCorner != myFS->SPH_VOLMAX) {
		oldIteration = ite;
		oldKCorr = kCorr;
		oldViscosity = visc;
		oldVorticity = vorticity;
		myFS->SPH_VOLMIN = currentMinCorner;
		myFS->SPH_VOLMAX = currentMaxCorner;
		myFS->setParameters(ite, visc, vorticity, kCorr);
		myFS->SPH_CreateExample(fluidPs);
		clearOrNot = true;
	}
	//runSimulation(currframe);

    UT_Interrupt *boss;
    if (error() < UT_ERROR_ABORT) {
		boss = UTgetInterrupt();
		gdp->clearAndDestroy();

		if (boss->opStart("Building Fluid") && 
			currframe < totalPos.size()) {	// currframe generation might not be able to catch up?
			for (auto& f : totalPos.at(currframe)) {
				/*UT_Vector3 pos;
				pos(0) = f.x;
				pos(1) = f.z;
				pos(2) = f.y;

				GU_PrimSphereParms sphere(gdp);
				double scale = myFS->SPH_RADIUS * 0.4;
				sphere.xform.scale(scale, scale, scale);
				sphere.xform.translate(pos);
				GU_PrimSphere::build(sphere, GEO_PRIMSPHERE);*/

				UT_Vector3 pos;
				pos(0) = f.x;
				pos(1) = f.z;
				pos(2) = f.y;

				GA_Offset ptoffstart = gdp->appendPoint();

				gdp->setPos3(ptoffstart, pos);
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

OP_ERROR SOP_Fluid::buildGeo() // this one not working
{
	UT_Interrupt* boss;
	if (error() < UT_ERROR_ABORT) {
		boss = UTgetInterrupt();
		gdp->clearAndDestroy();

		if (boss->opStart("Building Fluid") && currentFrame < totalPos.size()) {	// currframe generation might not be able to catch up?
			for (auto& f : totalPos.at(currentFrame)) {
				UT_Vector3 pos;
				pos(0) = f.x;
				pos(1) = f.z;
				pos(2) = f.y;

				GA_Offset ptoffstart = gdp->appendPoint();

				gdp->setPos3(ptoffstart, pos);
			}

			select(GU_SPrimitive);
		}

		boss->opEnd();
	}

	return error();
}