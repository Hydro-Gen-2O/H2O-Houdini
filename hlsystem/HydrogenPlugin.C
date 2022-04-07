


//#include <RE/RE_EGLServer.h>
#include <UT/UT_DSOVersion.h>
#include <SIM/SIM_DataFilter.h>
#include <SIM/SIM_Relationship.h>
#include <SIM/SIM_RelationshipGroup.h>
#include <OP/OP_OperatorTable.h>
#include <DOP/DOP_PRMShared.h>
#include <DOP/DOP_InOutInfo.h>
#include <DOP/DOP_Operator.h>
#include <DOP/DOP_Engine.h>

//LsystemPlugin
#include <UT/UT_Math.h>
#include <UT/UT_Interrupt.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <CH/CH_LocalVariable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_SpareData.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>


#include <limits.h>
#include "HydrogenPlugin.h"
using namespace HDK_Sample;

//
// Help is stored in a "wiki" style text file. 
//
// See the sample_install.sh file for an example.


///
/// newSopOperator is the hook that Houdini grabs from this dll
/// and invokes to register the SOP.  In this case we add ourselves
/// to the specified operator table.
///
void
newDopOperator(OP_OperatorTable *table)
{
	OP_Operator* op;

	op = new DOP_Operator("CusFluid",			// Internal name
		"Hydro-Gen2O",			// UI name
		DOP_Hydrogen::myConstructor,	// How to build the SOP
		DOP_Hydrogen::myTemplateList,	// My parameters
		1, //0,				// Min # of sources
		9999, //0,				// Max # of sources
		DOP_Hydrogen::myVariables,	// Local variables
		0,     // flag of operator
		1		//number of outputs (up to 4) from this operator
	);
	table->addOperator(op);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//PUT YOUR CODE HERE
//You need to declare your parameters here
//Example to declare a variable for angle you can do like this :
//static PRM_Name		angleName("angle", "Angle");
static PRM_Name		constraintIteration("constraintIteration", "Constraint Iteration");
static PRM_Name		constraintStiffness("constraintStiffness", "Constraint Stiffness");
static PRM_Name		viscosity("viscosity", "Viscosity");
static PRM_Name		vorticityConfinement("vorticityConfinement", "Vorticity Confinement");


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				     ^^^^^^^^    ^^^^^^^^^^^^^^^
//				     internal    descriptive version


// PUT YOUR CODE HERE
// You need to setup the initial/default values for your parameters here
// For example : If you are declaring the inital value for the angle parameter
//static PRM_Default angleDefault(30.0);	
static PRM_Default constraintIterationDefault(4);
static PRM_Default constraintStiffnessDefault(0.1);
static PRM_Default viscosityDefault(0.001);
static PRM_Default vorticityConfinementDefault(0.001);



static PRM_Name		 theInputIndexName("inputindex", "Input Index");
////////////////////////////////////////////////////////////////////////////////////////

PRM_Template
DOP_Hydrogen::myTemplateList[] = {
// PUT YOUR CODE HERE
// You now need to fill this template with your parameter name and their default value
// EXAMPLE : For the angle parameter this is how you should add into the template
//PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &angleName, &angleDefault, 0),
PRM_Template(PRM_INT,	PRM_Template::PRM_EXPORT_MIN, 1, &constraintIteration, &constraintIterationDefault, 0),
PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &constraintStiffness, &constraintStiffnessDefault, 0),
PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &viscosity, &viscosityDefault, 0),
PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &vorticityConfinement, &vorticityConfinementDefault, 0),


//code from DOP_GorupAndApply
// 
// Standard activation parameter.
   PRM_Template(PRM_INT_J,	1, &DOPactivationName,
			   &DOPactivationDefault),
	// Standard group parameter with group menu.
	PRM_Template(PRM_STRING,	1, &DOPgroupName, &DOPgroupDefault,
				&DOPgroupMenu),
	// The input index that determines which data will be attached to
	// each object.
	PRM_Template(PRM_INT_J,	1, &theInputIndexName, PRMzeroDefaults),

/////////////////////////////////////////////////////////////////////////////////////////////

    PRM_Template()
};


// Here's how we define local variables for the SOP.
enum {
	VAR_PT,		// Point number of the star
	VAR_NPT		// Number of points in the star
};

CH_LocalVariable
DOP_Hydrogen::myVariables[] = {
    { "PT",	VAR_PT, 0 },		// The table provides a mapping
    { "NPT",	VAR_NPT, 0 },		// from text string to integer token
    { 0, 0, 0 },
};

bool
DOP_Hydrogen::evalVariableValue(fpreal &val, int index, int thread)
{
    // myCurrPoint will be negative when we're not cooking so only try to
    // handle the local variables when we have a valid myCurrPoint index.
    if (myCurrPoint >= 0)
    {
	// Note that "gdp" may be null here, so we do the safe thing
	// and cache values we are interested in.
	switch (index)
	{
	    case VAR_PT:
		val = (fpreal) myCurrPoint;
		return true;
	    case VAR_NPT:
		val = (fpreal) myTotalPoints;
		return true;
	    default:
		/* do nothing */;
	}
    }
    // Not one of our variables, must delegate to the base class.
    return DOP_Node::evalVariableValue(val, index, thread);
}

OP_Node *
DOP_Hydrogen::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new DOP_Hydrogen(net, name, op);
}

DOP_Hydrogen::DOP_Hydrogen(OP_Network *net, const char *name, OP_Operator *op)
	: DOP_Node(net, name, op)
{
    myCurrPoint = -1;	// To prevent garbage values from being returned
}

DOP_Hydrogen::~DOP_Hydrogen() {}



//
// 
// 
// code from DOP_GroupAndApply.c
// 
// 
//

void
DOP_Hydrogen::processObjectsSubclass(fpreal time, int,
	const SIM_ObjectArray& objects,
	DOP_Engine& engine)
{
	SIM_ObjectArray	 filtered;
	UT_String		 group;
	int			 i, inputindex;

	// Grab the group string and filter our incoming objects using that
	// string. This narrows down the set of objects that we actually want
	// to operate on. The filtered array will contain only those objects
	// from the original objects array that match the supplied string.
	GROUP(group, time);
	SIM_DataFilterRootData	 filter(group);
	objects.filter(filter, filtered);

	// Loop through all the objects that passed the filter.
	for (i = 0; i < filtered.entries(); i++)
	{
		// Set information about the object we are going to process.
		// The first argument is the index of the current object within the
		// full list of objects we are going to process. The second
		// argument is the total number of objects we are going to process.
		// The last argument is a pointer to the actual object we are
		// processing.
		setCurrentObject(i, filtered.entries(), filtered(i));

		// The isActive function checks both the bypass flag and the
		// activation parameter on the node (if there is one, which there
		// is in this case). We call this function after calling
		// setCurrentObject and we call it for each object in case the
		// activation parameter uses some object-specific variables
		// like OBJID in an expression.
		if (isActive(time))
		{
			// Evaluate the input index. Also called after setCurrentObject
			// to properly evaluate objects specific local variables. Then
			// make sure there is something connected to the requested input.
			// We add one to the returned value to skip over the object
			// input.
			inputindex = INPUTINDEX(time) + 1;
			if (inputindex != 0 && getInput(inputindex))
			{
				SIM_Relationship* relationship;
				UT_String		 addtogroup;
				char			 numstr[UT_NUMBUF];

				// This function attaches the data connected to input
				// number inputindex to the current object.
				applyDataFromInput(time, inputindex, inputindex,
					*filtered(i),
					0, engine, 0, true);
				// Create a group name based on the input index.
				UT_String::itoa(numstr, inputindex);
				addtogroup = "applygroup_";
				addtogroup += numstr;

				// Create the relationship if it doesn't already exist,
				// and add the object to the group.
				relationship = engine.addRelationship(addtogroup,
					SIM_DATA_RETURN_EXISTING);
				if (relationship)
				{
					relationship->addGroup(filtered(i));
					SIM_DATA_CREATE(*relationship, SIM_RELGROUP_DATANAME,
						SIM_RelationshipGroup,
						SIM_DATA_RETURN_EXISTING);
				}
				else
					addError(DOP_CANTCREATERELATIONSHIP);
			}
		}
	}
}

void
DOP_Hydrogen::getInputInfoSubclass(int inputidx,
	DOP_InOutInfo& info) const
{
	// Our first input is an object input.
	// Our remaining inputs are data inputs.
	if (inputidx == 0)
		info = DOP_InOutInfo(DOP_INOUT_OBJECTS, false);
	else
		info = DOP_InOutInfo(DOP_INOUT_DATA, true);
}

void
DOP_Hydrogen::getOutputInfoSubclass(int outputidx,
	DOP_InOutInfo& info) const
{
	// Our single output is an object output.
	info = DOP_InOutInfo(DOP_INOUT_OBJECTS, false);
}

void
DOP_Hydrogen::GROUP(UT_String& str, fpreal t)
{
	evalString(str, DOPgroupName.getToken(), 0, t);
}

int
DOP_Hydrogen::INPUTINDEX(fpreal t)
{
	return evalInt(theInputIndexName.getToken(), 0, t);
}


//OP_ERROR
//SOP_Lsystem::cookMySop(OP_Context &context)
//{
//	fpreal		 now = context.getTime();
//
//	// PUT YOUR CODE HERE
//	// Decare the necessary variables and get always keep getting the current value in the node
//	// For example to always get the current angle thats set in the node ,you need to :
//	//    float angle;
//	//    angle = ANGLE(now)       
//    //    NOTE : ANGLE is a function that you need to use and it is declared in the header file to update your values instantly while cooking 
//	LSystem myplant;
//
//	float angle;
//	angle = ANGLE(now);
//
//	float stepSize;
//	stepSize = STEPSIZE(now);
//
//	int iterations;
//	iterations = ITERATIONS(now);
//
//	UT_String grammarFileUT;
//	evalString(grammarFileUT, "grammarFile", 0, now);
//	std::string grammarFile;
//	grammarFile = grammarFileUT.toStdString();
//
//	cout << "angle " << angle << endl;
//	cout << "stepSize " << stepSize << endl;
//	cout << "iterations " << iterations << endl;
//	cout << "grammar " << grammarFile << "\n" << endl;
// 
//
//	///////////////////////////////////////////////////////////////////////////
//
//	//PUT YOUR CODE HERE
//	// Next you need to call your Lystem cpp functions 
//	//Below is an example , you need to call the same functions based on the variables you declare
//     /*myplant.loadProgramFromString("F\nF->F[+F]F[-F]");  
//     myplant.setDefaultAngle(30.0f);
//     myplant.setDefaultStep(1.0f);*/
//	myplant.loadProgram(grammarFile);
//	myplant.setDefaultAngle(angle);
//	myplant.setDefaultStep(stepSize);
//
//
//
//	///////////////////////////////////////////////////////////////////////////////
//
//	// PUT YOUR CODE HERE
//	// You the need call the below function for all the genrations ,so that the end points points will be
//	// stored in the branches vector , you need to declare them first
//
//	
//	std::vector<LSystem::Branch> branches = std::vector<LSystem::Branch>();
//	myplant.process(iterations, branches);
//
//
//
//
//
//	///////////////////////////////////////////////////////////////////////////////////
//
//
//	// Now that you have all the branches ,which is the start and end point of each point ,its time to render 
//	// these branches into Houdini 
//    
//
//	// PUT YOUR CODE HERE
//	// Declare all the necessary variables for drawing cylinders for each branch 
//    float		 rad, tx, ty, tz;
//    int			 divisions, plane;
//    int			 xcoord =0, ycoord = 1, zcoord =2;
//    float		 tmp;
//    UT_Vector4		 pos;
//    GU_PrimPoly		*poly;
//    int			 i;
//    UT_Interrupt	*boss;
//
//
//    // Since we don't have inputs, we don't need to lock them.
//
//    divisions  = 5;	// We need twice our divisions of points
//    myTotalPoints = divisions;		// Set the NPT local variable value
//    myCurrPoint   = 0;			// Initialize the PT local variable
//
//    // Check to see that there hasn't been a critical error in cooking the SOP.
//    if (error() < UT_ERROR_ABORT)
//    {
//		boss = UTgetInterrupt();
//		if (divisions < 4)
//		{
//			// With the range restriction we have on the divisions, this
//			//	is actually impossible, but it shows how to add an error
//			//	message or warning to the SOP.
//			addWarning(SOP_MESSAGE, "Invalid divisions");
//			divisions = 4;
//		}
//		gdp->clearAndDestroy();
//
//		// Start the interrupt server
//		if (boss->opStart("Building LSYSTEM"))
//		{
//			// PUT YOUR CODE HERE
//			// Build a polygon
//			// You need to build your cylinders inside Houdini from here
//			// TIPS:
//			// Use GU_PrimPoly poly = GU_PrimPoly::build(see what values it can take)
//			// Also use GA_Offset ptoff = poly->getPointOffset()
//			// and gdp->setPos3(ptoff,YOUR_POSITION_VECTOR) to build geometry.
//			/*
//			UT_Vector3 fluidPos;
//			fluidPos(0) = myFS.points[0][0];
//			fluidPos(1) = myFS.points[0][1];
//			fluidPos(2) = myFS.points[0][2];
//
//			poly = GU_PrimPoly::build(gdp, 1, GU_POLY_CLOSED);
//
//			GA_Offset ptoffstart = poly->getPointOffset(0);
//			gdp->setPos3(ptoffstart, fluidPos);*/
//
//			for (int i = 0; i < branches.size(); i++)
//			{
//				UT_Vector3 posStart;
//				posStart(0) = branches.at(i).first[0];
//				posStart(1) = branches.at(i).first[2];
//				posStart(2) = branches.at(i).first[1];
//
//				UT_Vector3 posEnd;
//				posEnd(0) = branches.at(i).second[0];
//				posEnd(1) = branches.at(i).second[2];
//				posEnd(2) = branches.at(i).second[1];
//
//				cout << posStart[0] << " " << posStart[1] << " " << posStart[2] << endl;
//				cout << posEnd[0] << " " << posEnd[1] << " " << posEnd[2] << endl;
//
//				poly = GU_PrimPoly::build(gdp, 2, GU_POLY_CLOSED);
//
//				GA_Offset ptoffstart = poly->getPointOffset(0);
//				gdp->setPos3(ptoffstart, posStart);
//
//				GA_Offset ptoffend = poly->getPointOffset(1);
//				gdp->setPos3(ptoffend, posEnd);
//			}
//			////////////////////////////////////////////////////////////////////////////////////////////
//
//			// Highlight the star which we have just generated.  This routine
//			// call clears any currently highlighted geometry, and then it
//			// highlights every primitive for this SOP. 
//			select(GU_SPrimitive);
//		}
//
//		// Tell the interrupt server that we've completed. Must do this
//		// regardless of what opStart() returns.
//		boss->opEnd();
//    }
//
//    myCurrPoint = -1;
//    return error();
//}

