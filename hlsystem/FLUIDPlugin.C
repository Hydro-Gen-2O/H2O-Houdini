#include <UT/UT_DSOVersion.h>
//#include <RE/RE_EGLServer.h>

#include <UT/UT_Math.h>
#include <UT/UT_Interrupt.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <CH/CH_LocalVariable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_SpareData.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>

#include <GU/GU_PrimSphere.h>

#include <limits.h>
#include "FLUIDPlugin.h"

using namespace std;
using namespace HDK_Sample;

#define FRAME_RANGE 100

// newSopOperator is the hook that Houdini grabs from this dll
// and invokes to register the SOP.  In this case we add ourselves
// to the specified operator table.
void newSopOperator(OP_OperatorTable *table) {
    table->addOperator(
	    new OP_Operator("CusFluid",			// Internal name
			    "Hydro-Gen2O",			// UI name
			     SOP_Fluid::myConstructor,	// How to build the SOP
			     SOP_Fluid::myTemplateList,	// My parameters
			     0,				// Min # of sources
			     0,				// Max # of sources
			     SOP_Fluid::myVariables,	// Local variables
			     OP_FLAG_GENERATOR)		// Flag it as generator
	    );
}

//declare your parameters here
static PRM_Name		constraintIteration("constraintIteration", "Constraint Iteration");
static PRM_Name		artificialPressure("artificialPressure", "Artificial Pressure");
static PRM_Name		viscosity("viscosity", "Viscosity");
static PRM_Name		vorticityConfinement("vorticityConfinement", "Vorticity Confinement");
static PRM_Name		timeFrame("timeFrame", "Time Frame");
//				     ^^^^^^^^    ^^^^^^^^^^^^^^^
//				     internal    descriptive version

// initial/default values for your parameters here
static PRM_Default constraintIterationDefault(2);
static PRM_Default artificialPressureDefault(0.0001);
static PRM_Default viscosityDefault(0.01);
static PRM_Default vorticityConfinementDefault(0.0003);
static PRM_Default timeFrameDefault(0);

static PRM_Range iterationRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 30);
static PRM_Range tensileRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.01);
static PRM_Range viscosityRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.1);
static PRM_Range vorticityRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, 0.001);
static PRM_Range timeFrameRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_RESTRICTED, FRAME_RANGE);

PRM_Template
SOP_Fluid::myTemplateList[] = {
// default vals
PRM_Template(PRM_INT,	PRM_Template::PRM_EXPORT_MIN, 1, &constraintIteration, &constraintIterationDefault, 0, &iterationRange),
PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &artificialPressure, &artificialPressureDefault, 0, &tensileRange),
PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &viscosity, &viscosityDefault, 0, &viscosityRange),
PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &vorticityConfinement, &vorticityConfinementDefault, 0, &vorticityRange),
PRM_Template(PRM_INT,	PRM_Template::PRM_EXPORT_MIN, 1, &timeFrame, &timeFrameDefault, 0, &timeFrameRange),
    PRM_Template()
};

// Here's how we define local variables for the SOP.
enum {
	VAR_PT,		// Point number
	VAR_NPT,	// Number of points
	VAR_FS,  // my pointer to my fluid system
	VAR_OI, //iteration
	VAR_OK, //artificial pressure kCorr
	VAR_OVISC, //viscosity
	VAR_OVOR, //vorticity
};

CH_LocalVariable
SOP_Fluid::myVariables[] = {
    { "PT",	VAR_PT, 0 },		// The table provides a mapping
    { "NPT",	VAR_NPT, 0 },		// from text string to integer token
	{ "OI",	VAR_OI, 0 },
	{ "OK",	VAR_OK, 0 },
	{ "OVISC",	VAR_OVISC, 0 },
	{ "OVOR",	VAR_OVOR, 0 },
	{ "FS",	VAR_FS },
    { 0, 0, 0 },
};

bool SOP_Fluid::evalVariableValue(fpreal &val, int index, int thread) {
    // myCurrPoint will be negative when we're not cooking so only try to
    // handle the local variables when we have a valid myCurrPoint index.
    if (myCurrPoint >= 0) {
	// Note that "gdp" may be null here, so we do the safe thing
	// and cache values we are interested in.
	switch (index) {
	    case VAR_PT:
		val = (fpreal) myCurrPoint;
		return true;
	    case VAR_NPT:
		val = (fpreal) myTotalPoints;
		return true;
		case VAR_OI:
			val = (fpreal)oldIteration;
			return true;
		case VAR_OK:
			val = (fpreal)oldKCorr;
			return true;
		case VAR_OVISC:
			val = (fpreal)oldViscosity;
			return true;
		case VAR_OVOR:
			val = (fpreal)oldVorticity;
			return true;
	    default:
		/* do nothing */;
	}
    }
    // Not one of our variables, must delegate to the base class.
    return SOP_Node::evalVariableValue(val, index, thread);
}

OP_Node* SOP_Fluid::myConstructor(OP_Network *net, const char *name, OP_Operator *op) {
    return new SOP_Fluid(net, name, op);
}

SOP_Fluid::SOP_Fluid(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op) {
    myCurrPoint = -1;	// To prevent garbage values from being returned
	myFS = new FluidSystem();
	myFS->SPH_CreateExample(0, 0);

	//for test
	runSimulation(FRAME_RANGE);
	oldIteration = 2;
	oldKCorr = 0.0001;
	oldViscosity = 0.01;
	oldVorticity = 0.0003;
}

void SOP_Fluid::runSimulation(int frameNumber) {
	totalPos.clear();
	for (int i = 0; i <= frameNumber; ++i)
	{
		myFS->Run();

		std::vector<glm::dvec3> temp;
		for (auto& f : myFS->fluidPs) {
			glm::dvec3 scaledPos = f->pos;
			scaledPos /= SPH_RADIUS;
			temp.push_back(scaledPos);
		}
		totalPos.push_back(temp);
	}
}

SOP_Fluid::~SOP_Fluid() {}

unsigned SOP_Fluid::disableParms() {
    return 0;
}

OP_ERROR
SOP_Fluid::cookMySop(OP_Context &context)
{
	fpreal now = context.getTime();
	//FluidSystem myplant;

	int ite;
	ite = CONSTRAINT_ITERATION(now);

	float kCorr;
	kCorr = ARTIFICIAL_PRESSURE(now);

	float visc;
	visc = VISCOSITY(now);

	float vorticity;
	vorticity = VORTICITY_CONFINEMENT(now);

	int tf;
	tf = TIME_FRAME(now);

	if (ite != oldIteration || kCorr != oldKCorr || visc != oldViscosity || vorticity != oldVorticity)
	{
		oldIteration = ite;
		oldKCorr = kCorr;
		oldViscosity = visc;
		oldVorticity = vorticity;
		myFS->setParameters(ite, visc, vorticity, kCorr);

		myFS->SPH_CreateExample(0, 0);
		runSimulation(FRAME_RANGE);
	}

	/*std::cout << "iteration " << ite << endl;
	cout << "stiffnes " << stif << endl;
	cout << "viscosity " << visc << endl;
	cout << "vorticity " << vorticity << "\n" << endl;*/

	//PUT YOUR CODE HERE
	// Next you need to call your Lystem cpp functions 
	////Below is an example , you need to call the same functions based on the variables you declare
	//myplant.loadProgram(grammarFile);
	//myplant.setDefaultAngle(angle);
	//myplant.setDefaultStep(stepSize);

	// You the need call the below function for all the genrations ,so that the end points points will be
	// stored in the branches vector , you need to declare them first
	
	/*std::vector<LSystem::Branch> branches = std::vector<LSystem::Branch>();
	myplant.process(iterations, branches);*/

	//myFS->Run();
    float		 rad, tx, ty, tz;
    int			 divisions, plane;
    int			 xcoord =0, ycoord = 1, zcoord =2;
    float		 tmp;
    UT_Vector4	 pos;
    GU_PrimPoly		*poly;
    int			 i;
    UT_Interrupt	*boss;

    // Since we don't have inputs, we don't need to lock them.
    divisions = 5;	// We need twice our divisions of points
    myTotalPoints = divisions;		// Set the NPT local variable value
    myCurrPoint = 0;			// Initialize the PT local variable

    // Check to see that there hasn't been a critical error in cooking the SOP.
    if (error() < UT_ERROR_ABORT)
    {
		boss = UTgetInterrupt();
		//if (divisions < 4)
		//{
		//	// With the range restriction we have on the divisions, this
		//	//	is actually impossible, but it shows how to add an error
		//	//	message or warning to the SOP.
		//	addWarning(SOP_MESSAGE, "Invalid divisions");
		//	divisions = 4;
		//}
		gdp->clearAndDestroy();

		// Start the interrupt server
		if (boss->opStart("Building Fluid"))
		{
			// PUT YOUR CODE HERE
			// Build a polygon
			// You need to build your cylinders inside Houdini from here
			// TIPS:
			// Use GU_PrimPoly poly = GU_PrimPoly::build(see what values it can take)
			// Also use GA_Offset ptoff = poly->getPointOffset()
			// and gdp->setPos3(ptoff,YOUR_POSITION_VECTOR) to build geometry.

			for (auto& f : totalPos[tf]) {
				glm::dvec3 scaledPos = f;
				//scaledPos /= SPH_RADIUS;

				//debug
				//std::cout << "pos: " << scaledPos.x << " " << scaledPos.y << " " << scaledPos.z << std::endl;

				/*glTranslatef(scaledPos.x, scaledPos.y, scaledPos.z);
				glScalef(SPH_RADIUS, SPH_RADIUS, SPH_RADIUS);
				glColor4f(RED(f->clr), GRN(f->clr), BLUE(f->clr), ALPH(f->clr));
				drawSphere();*/


				UT_Vector3 pos;
				pos(0) = scaledPos.x;
				pos(1) = scaledPos.z;
				pos(2) = scaledPos.y;

				GU_PrimSphereParms sphere(gdp);
				sphere.xform.scale(SPH_RADIUS, SPH_RADIUS, SPH_RADIUS);
				sphere.xform.translate(pos);
				GU_PrimSphere::build(sphere, GEO_PRIMSPHERE);

				/*UT_Vector3 posUp;
				posUp(0) = scaledPos.x;
				posUp(1) = scaledPos.z + 0.2;
				posUp(2) = scaledPos.y;

				poly = GU_PrimPoly::build(gdp, 2, GU_POLY_CLOSED);
				GA_Offset ptoffstart = poly->getPointOffset(0);
				gdp->setPos3(ptoffstart, pos);

				GA_Offset ptoffend = poly->getPointOffset(1);
				gdp->setPos3(ptoffend, posUp);*/
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
    myCurrPoint = -1;
    return error();
}

