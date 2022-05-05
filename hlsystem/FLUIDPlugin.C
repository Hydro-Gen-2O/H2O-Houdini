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
static PRM_Name		PRM_viscosity("viscosity", "Viscosity");
static PRM_Name		vorticityConfinement("vorticityConfinement", "Vorticity Confinement");
static PRM_Name		PRM_minCorner("minCorner", "Min Corner");
static PRM_Name		PRM_maxCorner("maxCorner", "Max Corner");
//static PRM_Name		maxPts("maxPts", "Maximum Points");
static PRM_Name		simulateButton("simulateButton", "Run Simulation");
static PRM_Name		framesToBake("frameToBake", "Frames To Bake");
static PRM_Name		PRM_force("force", "Force");
//				     ^^^^^^^^    ^^^^^^^^^^^^^^^
//				     internal    descriptive version

// initial/default values for your parameters here
static PRM_Default constraintIterationDefault(2);
static PRM_Default artificialPressureDefault(0.0001);
static PRM_Default viscosityDefault(0.01);
static PRM_Default vorticityConfinementDefault(0.0000);
static PRM_Default frameBakeDefault(60);
static PRM_Default minDefault[] = { PRM_Default(-10.0), PRM_Default(0.0), PRM_Default(-10.0) };
static PRM_Default maxDefault[] = { PRM_Default(10.0), PRM_Default(20.0), PRM_Default(10.0) };
static PRM_Default forceDefault[] = { PRM_Default(0.0), PRM_Default(-9.8), PRM_Default(0.0) };
//static PRM_Default maxPtsDefault(5000);

static PRM_Range iterationRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 30);
static PRM_Range tensileRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.01);
static PRM_Range viscosityRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.1);
static PRM_Range vorticityRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.001);
static PRM_Range frameBakeRange(PRM_RANGE_RESTRICTED, 1, PRM_RANGE_RESTRICTED, 1000);
//static PRM_Range maxPtsRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 100000);

PRM_Template
SOP_Fluid::myTemplateList[] = {
	// default vals
	PRM_Template(PRM_INT,	PRM_Template::PRM_EXPORT_MIN, 1, &constraintIteration, &constraintIterationDefault, 0, &iterationRange),
	PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &artificialPressure, &artificialPressureDefault, 0, &tensileRange),
	PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &PRM_viscosity, &viscosityDefault, 0, &viscosityRange),
	PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &vorticityConfinement, &vorticityConfinementDefault, 0, &vorticityRange),
	PRM_Template(PRM_XYZ_J, 3, &PRM_minCorner, minDefault),
	PRM_Template(PRM_XYZ_J, 3, &PRM_maxCorner, maxDefault),
	PRM_Template(PRM_XYZ_J, 3, &PRM_force, forceDefault),
	PRM_Template(PRM_INT,	PRM_Template::PRM_EXPORT_MIN, 1, &framesToBake, &frameBakeDefault, 0, &frameBakeRange),
	//PRM_Template(PRM_INT,	PRM_Template::PRM_EXPORT_MIN, 1, &maxPts, &maxPtsDefault, 0, &maxPtsRange),
	PRM_Template(PRM_CALLBACK, 1, &simulateButton, 0, 0, 0, &simulate),
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
	validFluidPs = true;
}

int SOP_Fluid::simulate(void* op, int index, fpreal t, const PRM_Template*) {
	SOP_Fluid* fluid = (SOP_Fluid*)op;
	if (fluid->validFluidPs) {
		fluid->runSimulation(fluid->currentFrame, true);
		return 1;
	}
	return -1;
}

void SOP_Fluid::runSimulation(int frameNumber, bool refresh) {
	if (refresh) {
		totalPos.clear();
		myFS->setParameters(iters, viscosity, vorticity, kcorr);
		myFS->SPH_VOLMIN = minCorner;
		myFS->SPH_VOLMAX = maxCorner;
		myFS->FORCE = force;
		myFS->SPH_CreateExample(fluidPs);
	}
	int moreFrame = frameNumber - totalPos.size();
	if (moreFrame >= 0) {
		for (int i = 0; i <= moreFrame + frameRange; ++i) {
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

	if (!validFluidPs) {
		addWarning(SOP_MESSAGE, "Fluid volume out of bounds! Decrease fluid volume or increase bounds.");
	}

	// 	// update the simulation values - but do not run, only run on callback - from button
	iters = CONSTRAINT_ITERATION(now);
	kcorr = ARTIFICIAL_PRESSURE(now);
	viscosity = VISCOSITY(now);
	vorticity = VORTICITY_CONFINEMENT(now);
	float minx = evalFloat("minCorner", 0, now);
	float miny = evalFloat("minCorner", 1, now);
	float minz = evalFloat("minCorner", 2, now);
	float maxx = evalFloat("maxCorner", 0, now);
	float maxy = evalFloat("maxCorner", 1, now);
	float maxz = evalFloat("maxCorner", 2, now);
	float forcex = evalFloat("force", 0, now);
	float forcey = evalFloat("force", 1, now);
	float forcez = evalFloat("force", 2, now);
	force = glm::dvec3(forcex, forcez, forcey);
	frameRange = FRAME_BAKE(now);
	minCorner = glm::dvec3(minx, minz, miny); // flip z & y
	maxCorner = glm::dvec3(maxx, maxz, maxy);
	//int maxPts = MAX_PTS(now);

	//get inputs
	OP_AutoLockInputs inputs(this);
	if (inputs.lockInput(0, context) >= UT_ERROR_ABORT) { return error(); }
	// only if input geo is different, re-get all the pts. again, do not run - only callback runs
	int input_changed;
	duplicateChangedSource(0, context, &input_changed);
	if (input_changed) {
		fluidPs.clear();
		GU_Detail* fluid_gdp = new GU_Detail(inputGeo(0, context));
		GA_Offset ptoff;
		GA_FOR_ALL_PTOFF(fluid_gdp, ptoff) {
			UT_Vector3 pos = fluid_gdp->getPos3(ptoff);
			glm::dvec3 p(pos[0], pos[2], pos[1]);
			if (p.x > minCorner.x && p.y > minCorner.y && p.z > minCorner.z &&
				p.x < maxCorner.x && p.y < maxCorner.y && p.z < maxCorner.z) {
				fluidPs.push_back(p);
			} else {
				addWarning(SOP_MESSAGE, "Fluid volume out of bounds! Decrease fluid volume or increase bounds.");
				validFluidPs = false;
				return error();
			}
		}
		validFluidPs = true;
	}
	
	runSimulation(currframe, false); // update if user is scrubbing

    UT_Interrupt *boss;
    if (error() < UT_ERROR_ABORT) {
		boss = UTgetInterrupt();
		gdp->clearAndDestroy();

		if (boss->opStart("Building Fluid") && 
			currframe < totalPos.size()) {	// currframe generation might not be able to catch up
			for (auto& f : totalPos.at(currframe)) {
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