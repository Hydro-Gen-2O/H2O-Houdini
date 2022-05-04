#ifndef __FLUID_PLUGIN_h__
#define __FLUID_PLUGIN_h__

//#include <GEO/GEO_Point.h>
#include <SOP/SOP_Node.h>
#include "fluid_system.h"

class SOP_Fluid : public SOP_Node {
public:
    static OP_Node *myConstructor(OP_Network*, const char *name, OP_Operator *op);
    /// Stores the description of the interface of the SOP in Houdini.
    /// Each parm template refers to a parameter.
    static PRM_Template		 myTemplateList[];

    /// This optional data stores the list of local variables.
    static CH_LocalVariable	 myVariables[];
protected:
    SOP_Fluid(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_Fluid();

    /// cookMySop does the actual work of the SOP computing
    virtual OP_ERROR cookMySop(OP_Context &context);

    void runSimulation(int frameNumber, bool reRun);

    // callback used by the "Clear All" parameter
    static int simulate(void* op, int index, fpreal time, const PRM_Template*);
    OP_ERROR buildGeo();
private:
	// functions to constantly update the cook function, get the current value that the node has
    //exint MAX_PTS(exint t) { return evalFloat("maxPts", 0, t); }

    exint CONSTRAINT_ITERATION(exint t) { return evalInt("constraintIteration", 0, t); }
    exint FRAME_BAKE(exint t) { return evalInt("frameToBake", 0, t); }
    fpreal ARTIFICIAL_PRESSURE(fpreal t) { return evalFloat("artificialPressure", 0, t); }
    fpreal VISCOSITY(fpreal t) { return evalFloat("viscosity", 0, t); }
    fpreal VORTICITY_CONFINEMENT(fpreal t) { return evalFloat("vorticityConfinement", 0, t); }
    //exint START_FRAME(exint t) { return evalInt("startFrame", 0, t); }

    glm::dvec3 force;
    bool validFluidPs;
    int frameRange;
    int iters;
    //int     myStartFrame;
    float kcorr;
    float viscosity;
    float vorticity;
    glm::dvec3 maxCorner;
    glm::dvec3 minCorner;
    int     currentFrame; //for calculation in simulate

    OP_Context* myContext;
    FluidSystem* myFS;
    std::vector<std::vector<glm::dvec3>> totalPos; 
    std::vector<glm::dvec3> fluidPs;
};
#endif
