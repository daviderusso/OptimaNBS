import com.gurobi.gurobi.*;

import java.util.HashMap;
import java.util.ArrayList;
import java.util.Objects;

/**
 * SolverMILPGurobi is a class for solving Mixed-Integer Linear Programming (MILP) problems
 * using the Gurobi optimization library. It is designed for optimization problems related to
 * urban planning, including fairness, resource allocation, and clustering.
 */
public class SolverMILPGurobi {
    // Instance of the problem containing data and parameters
    Instance instance;

    // A set of kernels used for spatial computations, organized by type and name
    HashMap<String, HashMap<String, Kernel>> kernelSet;

    // Gurobi environment and model for optimization
    GRBEnv env;
    GRBModel model;

    // Time limit for the optimization process
    double timeLimit;

    // Reference values for the objectives
    double[] refVal;

    // Maximum possible values for the objectives
    double[] maxVal;

    // Bounds for the budget changes (delta)
    double[] deltaB;

    // Fairness bounds for the objectives
    double[] boundFairness;

    // A reference matrix used for certain constraints in the model
    double[][][] refMatrix;

    // A clustered map of parks, used for spatial constraints
    int[][] clusteredMapPark;

    // Size of each cluster in the clustered map
    ArrayList<Integer> clusterSize;

    // Configuration object containing model-specific settings
    Config config;

    /**
     * Constructor for SolverMILPGurobi.
     * Initializes the solver with the problem instance, kernel data, and optimization parameters.
     *
     * @param inst             the instance of the problem
     * @param kSet             the set of kernels used for calculations
     * @param showOut          a flag to enable or disable Gurobi output
     * @param timeLimit        the time limit for the optimization process
     * @param refVal           reference values for the objectives
     * @param deltaB           bounds for the budget changes
     * @param maxVal           maximum possible values for the objectives
     * @param clusteredMapPark clustered map of parks for spatial constraints
     * @param clusterSize      sizes of the clusters in the clustered map
     * @param boundFairness    fairness bounds for the objectives
     * @param config           configuration object with additional settings
     * @throws GRBException if an error occurs during the initialization of the Gurobi model or environment
     */
    public SolverMILPGurobi(
            Instance inst,
            HashMap<String, HashMap<String, Kernel>> kSet,
            boolean showOut,
            double timeLimit,
            double[] refVal,
            double[] deltaB,
            double[] maxVal,
            int[][] clusteredMapPark,
            ArrayList<Integer> clusterSize,
            double[] boundFairness,
            Config config
    ) throws GRBException {
        // Initialize instance, configuration, and kernel set
        this.instance = inst;
        this.config = config;
        this.kernelSet = kSet;

        // Create a Gurobi environment and model
        this.env = new GRBEnv();
        this.model = new GRBModel(env);

        // Set the Gurobi output flag based on the showOut parameter
        if (!showOut) {
            model.set(GRB.IntParam.OutputFlag, 0);
        }

        // Set the time limit for optimization
        this.timeLimit = timeLimit;

        // Initialize budget change bounds and fairness bounds
        this.deltaB = deltaB;
        this.boundFairness = boundFairness;

        // Initialize maximum and reference values for objectives
        this.maxVal = maxVal;
        this.refVal = refVal;

        // Initialize the reference matrix for each urban challenge
        this.refMatrix = new double[inst.UrbanChallenges.length][inst.W][inst.H];
        for (int b = 0; b < inst.UrbanChallenges.length; b++) {
            for (int i = 0; i < inst.W; i++) {
                for (int j = 0; j < inst.H; j++) {
                    refMatrix[b][i][j] = refVal[b]; // Populate with reference values
                }
            }
        }

        // Initialize spatial constraints using clustered maps
        this.clusteredMapPark = clusteredMapPark;
        this.clusterSize = clusterSize;
    }

    /**
     * Computes the Z value for a specific cell (i, j) based on the decision variables (xVar)
     * and the kernel functions (K_t). Excludes pre-existing installations from the computation.
     *
     * @param xVar a 3D array of decision variables representing potential installations
     * @param i    the row index of the cell
     * @param j    the column index of the cell
     * @param K_t  a mapping of kernel functions for each green type
     * @return a linear expression representing the Z value for the cell
     */
    private GRBLinExpr computeZ(GRBVar[][][] xVar, int i, int j, HashMap<String, Kernel> K_t) {
        GRBLinExpr val = new GRBLinExpr(); // Initialize the linear expression for Z
        double[][] K; // The kernel matrix
        int a; // Vertical kernel offset (half the kernel height)
        int b; // Horizontal kernel offset (half the kernel width)
        int W = xVar[0].length; // Grid width
        int H = xVar[0][0].length; // Grid height

        // Iterate over each green type
        for (int t = 0; t < instance.GreenType.length; t++) {
            K = K_t.get(instance.GreenType[t]).k; // Get the kernel for the current green type
            a = K.length / 2;
            b = K[0].length / 2;

            // Process the top-left quadrant of the kernel
            for (int ii = 0; ii < a + 1; ii++) {
                if (i - (a - ii) >= 0) { // Ensure the kernel stays within bounds
                    for (int jj = 0; jj < b + 1; jj++) {
                        if (j - (b - jj) >= 0) {
                            if (!instance.checkPreExist(t, i - (a - ii), j - (b - jj))) { // Exclude pre-existing installations
                                val.addTerm(K[ii][jj], xVar[t][i - (a - ii)][j - (b - jj)]);
                            }
                        }
                    }
                    // Process the top-right quadrant of the kernel
                    for (int jj = b + 1; jj < K.length; jj++) {
                        if (j + (jj - b) < H) {
                            if (!instance.checkPreExist(t, i - (a - ii), j + (jj - b))) { // Exclude pre-existing installations
                                val.addTerm(K[ii][jj], xVar[t][i - (a - ii)][j + (jj - b)]);
                            }
                        }
                    }
                }
            }

            // Process the bottom-left quadrant of the kernel
            for (int ii = a + 1; ii < K.length; ii++) {
                if (i + (ii - a) < W) {
                    for (int jj = 0; jj < b + 1; jj++) {
                        if (j - (b - jj) >= 0) {
                            if (!instance.checkPreExist(t, i + (ii - a), j - (b - jj))) { // Exclude pre-existing installations
                                val.addTerm(K[ii][jj], xVar[t][i + (ii - a)][j - (b - jj)]);
                            }
                        }
                    }
                    // Process the bottom-right quadrant of the kernel
                    for (int jj = b + 1; jj < K[0].length; jj++) {
                        if (j + (jj - b) < H) {
                            if (!instance.checkPreExist(t, i + (ii - a), j + (jj - b))) { // Exclude pre-existing installations
                                val.addTerm(K[ii][jj], xVar[t][i + (ii - a)][j + (jj - b)]);
                            }
                        }
                    }
                }
            }
        }

        return val; // Return the computed Z value as a linear expression
    }

    /**
     * Computes the Z value for a specific cell (i, j) considering fairness constraints.
     * Does not exclude pre-existing installations from the computation.
     *
     * @param xVar a 3D array of decision variables representing potential installations
     * @param i    the row index of the cell
     * @param j    the column index of the cell
     * @param K_t  a mapping of kernel functions for each green type
     * @return a linear expression representing the Z value for the cell with fairness constraints
     */
    private GRBLinExpr computeZFairness(GRBVar[][][] xVar, int i, int j, HashMap<String, Kernel> K_t) {
        GRBLinExpr val = new GRBLinExpr(); // Initialize the linear expression for Z
        double[][] K; // The kernel matrix
        int a; // Vertical kernel offset (half the kernel height)
        int b; // Horizontal kernel offset (half the kernel width)
        int W = xVar[0].length; // Grid width
        int H = xVar[0][0].length; // Grid height

        // Iterate over each green type
        for (int t = 0; t < instance.GreenType.length; t++) {
            K = K_t.get(instance.GreenType[t]).k; // Get the kernel for the current green type
            a = K.length / 2;
            b = K[0].length / 2;

            // Process the top-left quadrant of the kernel
            for (int ii = 0; ii < a + 1; ii++) {
                if (i - (a - ii) >= 0) { // Ensure the kernel stays within bounds
                    for (int jj = 0; jj < b + 1; jj++) {
                        if (j - (b - jj) >= 0) {
                            val.addTerm(K[ii][jj], xVar[t][i - (a - ii)][j - (b - jj)]);
                        }
                    }
                    // Process the top-right quadrant of the kernel
                    for (int jj = b + 1; jj < K.length; jj++) {
                        if (j + (jj - b) < H) {
                            val.addTerm(K[ii][jj], xVar[t][i - (a - ii)][j + (jj - b)]);
                        }
                    }
                }
            }

            // Process the bottom-left quadrant of the kernel
            for (int ii = a + 1; ii < K.length; ii++) {
                if (i + (ii - a) < W) {
                    for (int jj = 0; jj < b + 1; jj++) {
                        if (j - (b - jj) >= 0) {
                            val.addTerm(K[ii][jj], xVar[t][i + (ii - a)][j - (b - jj)]);
                        }
                    }
                    // Process the bottom-right quadrant of the kernel
                    for (int jj = b + 1; jj < K[0].length; jj++) {
                        if (j + (jj - b) < H) {
                            val.addTerm(K[ii][jj], xVar[t][i + (ii - a)][j + (jj - b)]);
                        }
                    }
                }
            }
        }
        return val; // Return the computed Z value as a linear expression
    }

    /**
     * This method computes the Mixed-Integer Linear Programming (MILP) model for optimizing the placement of green types in a grid
     * considering urban challenges, costs, and fairness. It initializes the decision variables, the objective function,
     * and the necessary constraints. The model is then set up to minimize the overall cost while optimizing for urban challenges
     * and maximizing fairness, subject to budget and placement constraints.
     * <p>
     * The decision variables include binary variables for placing green types and continuous variables for urban challenge metrics,
     * with auxiliary variables to linearize the model. Constraints ensure that only one green type is placed in each cell and that
     * the total cost does not exceed the available budget.
     *
     * @throws GRBException if there is an error during model creation or optimization with the Gurobi solver
     */
    public void compute_MILP() throws GRBException {
        // Declare a 3D array of GRBVar for decision variables representing the placement of different green types
        GRBVar[][][] x_var = new GRBVar[instance.GreenType.length][instance.W][instance.H];

        // Initialize binary variables x_var[t][i][j] indicating whether green type 't' is placed at position (i, j)
        for (int t = 0, l = instance.GreenType.length; t < l; t++) {
            for (int i = 0; i < instance.W; i++) {
                for (int j = 0; j < instance.H; j++) {
                    x_var[t][i][j] = this.model.addVar(0, 1, 0, GRB.BINARY, "x_" + t + "_" + i + "_" + j);
                }
            }
        }

        // Declare a 3D array of GRBVar for continuous variables z_var for results from applying the kernel
        GRBVar[][][] z_var = new GRBVar[instance.UrbanChallenges.length][instance.W][instance.H];

        // Initialize continuous variables z_var[b][i][j] for each urban challenge 'b' at position (i, j)
        for (int b = 0, l2 = instance.UrbanChallenges.length; b < l2; b++) {
            for (int i = 0; i < instance.W; i++) {
                for (int j = 0; j < instance.H; j++) {
                    z_var[b][i][j] = this.model.addVar(0.0, config.BigM[b], 0, GRB.CONTINUOUS, "z_" + b + "_" + i + "_" + j);
                }
            }
        }

        // Declare and initialize max value variables for each urban challenge
        GRBVar[] z_max = new GRBVar[instance.UrbanChallenges.length];
        for (int b = 0, l2 = instance.UrbanChallenges.length; b < l2; b++) {
            z_max[b] = this.model.addVar(0.0, 1.0, 0, GRB.CONTINUOUS, "z_max_" + b);
        }

        // Declare and initialize mean value variables for each urban challenge
        GRBVar[] z_mean = new GRBVar[instance.UrbanChallenges.length];
        for (int b = 0, l2 = instance.UrbanChallenges.length; b < l2; b++) {
            z_mean[b] = this.model.addVar(0.0, 1.0, 0, GRB.CONTINUOUS, "z_mean_" + b);
        }

        // Declare and initialize z_bar_var for kernel results in range for each urban challenge
        GRBVar[][][] z_bar_var = new GRBVar[instance.UrbanChallenges.length][instance.W][instance.H];
        for (int b = 0, l2 = instance.UrbanChallenges.length; b < l2; b++) {
            for (int i = 0; i < instance.W; i++) {
                for (int j = 0; j < instance.H; j++) {
                    z_bar_var[b][i][j] = this.model.addVar(0.0, deltaB[b], 0, GRB.CONTINUOUS, "z_bar_" + b + "_" + i + "_" + j);
                }
            }
        }

        // Declare and initialize auxiliary variables for linearizing the product of z and y
        GRBVar[][][] z_mu = new GRBVar[instance.UrbanChallenges.length][instance.W][instance.H];
        for (int b = 0, l2 = instance.UrbanChallenges.length; b < l2; b++) {
            for (int i = 0; i < instance.W; i++) {
                for (int j = 0; j < instance.H; j++) {
                    z_mu[b][i][j] = this.model.addVar(0.0, config.BigM[b], 0, GRB.CONTINUOUS, "z_mu_" + b + "_" + i + "_" + j);
                }
            }
        }

        // Declare a binary variable y_var to represent whether the reduction range is exceeded for each urban challenge
        GRBVar[][][] y_var = new GRBVar[instance.UrbanChallenges.length][instance.W][instance.H];
        for (int b = 0, l = instance.UrbanChallenges.length; b < l; b++) {
            for (int i = 0; i < instance.W; i++) {
                for (int j = 0; j < instance.H; j++) {
                    y_var[b][i][j] = this.model.addVar(0.0, 1.0, 0, GRB.BINARY, "y_" + b + "_" + i + "_" + j);
                }
            }
        }

        // Declare a continuous variable for fairness in the solution
        GRBVar fairness = this.model.addVar(0.0, config.BigM[config.BigM.length - 1], 0, GRB.CONTINUOUS, "fairness");


        ////////////////////////////////////////FO: Objective Function
        GRBLinExpr FO = new GRBLinExpr();
        int IdBmax = 0;
        int IdBavg = 0;

        // Minimize the urban challenge metrics (e.g., max or avg values)
        for (int b = 0, l2 = config.UCsOF.size() - 2; b < l2; b++) {
            if (config.UCsOF.get(b).contains(Config.Avg)) {
                FO.addTerm(Config.theta[b], z_mean[IdBavg]);
                IdBavg++;
            } else {
                FO.addTerm(Config.theta[b], z_max[IdBmax]);
                IdBmax++;
            }
        }

        // Minimize the costs associated with each green type in the grid (excluding pre-existing cells)
        GRBLinExpr obj = new GRBLinExpr();
        for (int t = 0, l = instance.GreenType.length; t < l; t++) {
            for (int i = 0; i < instance.W; i++) {
                for (int j = 0; j < instance.H; j++) {
                    if (!instance.checkPreExist(t, i, j)) {  // Exclude pre-existing variables
                        obj.addTerm(config.normalizedCost.get(instance.GreenType[t]), x_var[t][i][j]);
                    }
                }
            }
        }
        // Add the cost minimization term to the objective function
        FO.multAdd(Config.theta[Config.theta.length - 2], obj);

        // Maximize fairness by minimizing the opposite of fairness
        FO.addTerm((-1) * (Config.theta[Config.theta.length - 1]), fairness);

        // Set the objective to minimize the defined objective function (FO)
        model.setObjective(FO, GRB.MINIMIZE);

        ////////////////////////////////////////CONSTRAINTS:
        // Add constraints to ensure that only one green type is assigned per cell
        for (int i = 0; i < instance.W; i++) {
            for (int j = 0; j < instance.H; j++) {
                GRBLinExpr constrOneGreen = new GRBLinExpr();
                for (int t = 0, l = instance.GreenType.length; t < l; t++) {
                    constrOneGreen.addTerm(1, x_var[t][i][j]);
                }
                // Ensure that no more than one green type is placed in each cell
                model.addConstr(constrOneGreen, GRB.LESS_EQUAL, 1, "OneGreen_" + i + "_" + j);
            }
        }

        // Add budget constraint to ensure that the total cost does not exceed the budget
        GRBLinExpr constrBudget = new GRBLinExpr();
        for (int i = 0; i < instance.W; i++) {
            for (int j = 0; j < instance.H; j++) {
                for (int t = 0, l = instance.GreenType.length; t < l; t++) {
                    if (!instance.checkPreExist(t, i, j)) {  // Exclude pre-existing variables
                        constrBudget.addTerm(config.normalizedCost.get(instance.GreenType[t]), x_var[t][i][j]);
                    }
                }
            }
        }
        // Ensure that the total cost does not exceed the normalized budget
        model.addConstr(constrBudget, GRB.LESS_EQUAL, config.NormalizedBudget, "Budget");


        // Fixes pre-existing variables to 1 for each green type (e.g., if a green type is already placed at a location, it should remain fixed to that placement)
        for (int t = 0, len = instance.GreenType.length; t < len; t++) {
            for (int p = 0, lenP = instance.PreExistent.get(instance.GreenType[t]).size(); p < lenP; p++) {
                model.addConstr(x_var[t][instance.PreExistent.get(instance.GreenType[t]).get(p).r][instance.PreExistent.get(instance.GreenType[t]).get(p).c],
                        GRB.EQUAL, 1, "pre_exist_" + instance.PreExistent.get(instance.GreenType[t]).get(p).r + "-" + instance.PreExistent.get(instance.GreenType[t]).get(p).c);
            }
        }

        // Fixes forbidden variables to 0 for each green type (e.g., certain positions may be restricted for specific green types)
        int[][] mat_of_forbidden = new int[instance.W][instance.H];
        int curr_i;
        int curr_j;
        for (int t = 0, len = instance.GreenType.length; t < len; t++) {
            if (Objects.equals(instance.GreenType[t], Config.ClusteredGT)) {
                for (int i = 0; i < instance.W; i++) {
                    for (int j = 0; j < instance.H; j++) {
                        // Apply forbidden constraints based on the clusteredMapPark (0 means forbidden)
                        if (clusteredMapPark[i][j] == 0) {
                            model.addConstr(x_var[t][i][j], GRB.EQUAL, 0, "Forbidden_" + i + "_" + j + "_" + t);
                            mat_of_forbidden[i][j]++;
                        }
                    }
                }
            } else {
                // For non-clustered green types, apply specific forbidden constraints
                for (int f = 0, lenF = instance.Forbidden.get(instance.GreenType[t]).size(); f < lenF; f++) {
                    curr_i = instance.Forbidden.get(instance.GreenType[t]).get(f).r;
                    curr_j = instance.Forbidden.get(instance.GreenType[t]).get(f).c;
                    model.addConstr(x_var[t][curr_i][curr_j], GRB.EQUAL, 0, "Forbidden_" + curr_i + "_" + curr_j + "_" + t);
                    mat_of_forbidden[curr_i][curr_j]++;
                }
            }
        }

        // Kernel application constraints: Ensure that the z_var values for urban challenges are correctly calculated
        for (int b = 0, l2 = instance.UrbanChallenges.length; b < l2; b++) {
            for (int i = 0; i < instance.W; i++) {
                for (int j = 0; j < instance.H; j++) {
                    model.addConstr(z_var[b][i][j], GRB.EQUAL, this.computeZ(x_var, i, j, kernelSet.get(instance.UrbanChallenges[b])), "Kernel" + b + "_" + i + "_" + j);
                }
            }
        }

        // z_max >= A - z_bar constraint and calculation of z_mean as the average of (A - z_bar)
        for (int b = 0, l2 = instance.UrbanChallenges.length; b < l2; b++) {
            GRBLinExpr valCurr = new GRBLinExpr();
            for (int i = 0; i < instance.W; i++) {
                for (int j = 0; j < instance.H; j++) {
                    if ((mat_of_forbidden[i][j] < instance.GreenType.length) || (InstanceUtils.checkAdjTile(instance, mat_of_forbidden, i, j, kernelSet))) {
                        GRBLinExpr expr = new GRBLinExpr();
                        expr.addConstant(instance.A[b][i][j]);  // A is an array of challenge-specific values
                        expr.addTerm(-1.0, z_bar_var[b][i][j]); // Apply z_bar variable
                        model.addConstr(z_max[b], GRB.GREATER_EQUAL, expr, "Z_max" + b + "_" + i + "_" + j);
                        valCurr.multAdd(1.0, expr); // Accumulate the terms for z_mean calculation
                    }
                }
            }
            // z_mean is the average of (A - z_bar)
            GRBLinExpr valMean = new GRBLinExpr();
            valMean.multAdd((1.0 / (instance.W * instance.H)), valCurr);
            model.addConstr(valMean, GRB.EQUAL, z_mean[b], "Z_mean" + b);
        }

        // Set value for z_bar variable, linearizing the computation of kernel values
        for (int b = 0, l2 = instance.UrbanChallenges.length; b < l2; b++) {
            for (int i = 0; i < instance.W; i++) {
                for (int j = 0; j < instance.H; j++) {
                    // Ensure that z_bar is the difference between z_var and z_mu
                    GRBLinExpr expr = new GRBLinExpr();
                    expr.addTerm(1.0, z_var[b][i][j]);
                    expr.addTerm(-1.0, z_mu[b][i][j]);
                    model.addConstr(z_bar_var[b][i][j], GRB.EQUAL, expr, "z_bar_" + b + "_" + i + "_" + j);

                    // Constrain z_mu with respect to the BigM values and y_var (used for linearization)
                    GRBLinExpr expr1 = new GRBLinExpr();
                    expr1.addTerm(config.BigM[b], y_var[b][i][j]);
                    model.addConstr(z_mu[b][i][j], GRB.LESS_EQUAL, expr1, "z_mu" + b + "_" + i + "_" + j + "_LE");

                    GRBLinExpr expr2 = new GRBLinExpr();
                    expr2.addConstant(1.0);
                    expr2.addTerm(-1.0, y_var[b][i][j]);

                    GRBLinExpr expr3 = new GRBLinExpr();
                    expr3.multAdd((-1.0) * config.BigM[b], expr2);
                    expr3.addConstant((-1.0) * this.deltaB[b]);
                    expr3.addTerm(1.0, z_var[b][i][j]);
                    model.addConstr(z_mu[b][i][j], GRB.GREATER_EQUAL, expr3, "z_mu" + b + "_" + i + "_" + j + "_GE");
                }
            }
        }

        // Linking constraints between y and z variables
        for (int b = 0, l2 = instance.UrbanChallenges.length; b < l2; b++) {
            for (int i = 0; i < instance.W; i++) {
                for (int j = 0; j < instance.H; j++) {
                    // Define the constraints that link y_var and z_var to ensure consistency in the model
                    GRBLinExpr expr1 = new GRBLinExpr();
                    expr1.addTerm(1.0, z_var[b][i][j]);
                    expr1.addConstant((-1) * deltaB[b]);
                    GRBLinExpr expr2 = new GRBLinExpr();
                    expr2.addTerm(config.BigM[b], y_var[b][i][j]);
                    model.addConstr(expr1, GRB.LESS_EQUAL, expr2, "link1_z_var_" + b + "_" + i + "_" + j);

                    GRBLinExpr expr3 = new GRBLinExpr();
                    expr3.addConstant(deltaB[b]);
                    expr3.addTerm(-1.0, z_var[b][i][j]);
                    GRBLinExpr expr4 = new GRBLinExpr();
                    expr4.addConstant(1);
                    expr4.addTerm(-1, y_var[b][i][j]);
                    GRBLinExpr expr5 = new GRBLinExpr();
                    expr5.multAdd(config.BigM[b], expr4);
                    model.addConstr(expr3, GRB.LESS_EQUAL, expr5, "link2_z_var_" + b + "_" + i + "_" + j);
                }
            }
        }

        // Handle urban park cluster sizes
        int UPid = -1;
        for (int i = 0; i < instance.GreenType.length; i++) {
            if (Objects.equals(instance.GreenType[i], Config.ClusteredGT)) {
                UPid = i;
                break;
            }
        }

        // Define a binary variable to represent whether a cluster is used
        GRBVar[] clusters = new GRBVar[this.clusterSize.size()];
        for (int b = 0, l2 = this.clusterSize.size(); b < l2; b++) {
            if (clusterSize.get(b) > 0) {
                clusters[b] = this.model.addVar(0, 1, 0, GRB.BINARY, "cluster_" + b);
            }
        }

        // If a cluster is used, all locations within the cluster must be assigned the green type
        for (int c = 0; c < clusters.length; c++) {
            if (clusterSize.get(c) > 0) {
                for (int i = 0; i < instance.W; i++) {
                    for (int j = 0; j < instance.H; j++) {
                        if (clusteredMapPark[i][j] == c) {
                            model.addConstr(x_var[UPid][i][j], GRB.EQUAL, clusters[c], "cluster_" + c + "_" +  UPid + "_" + i + "_" + j);
                        }
                    }
                }
            }
        }

        // Calculate fairness based on population distribution and urban challenges
        GRBLinExpr fairnessSum = new GRBLinExpr();
        int idFair = instance.UrbanChallenges.length - 1;
        for (int i = 0; i < instance.W; i++) {
            for (int j = 0; j < instance.H; j++) {
                double currPop = this.instance.A[idFair][i][j]; // Population data for fairness calculation
                GRBLinExpr fairK = this.computeZFairness(x_var, i, j, kernelSet.get(instance.UrbanChallenges[idFair]));
                fairnessSum.multAdd(currPop, fairK); // Weighted sum based on population
            }
        }

        // Final fairness constraint ensuring the fairness value lies between the specified bounds
        fairnessSum.addConstant(-1.0 * boundFairness[0]);
        GRBLinExpr fairnessFinal = new GRBLinExpr();
        fairnessFinal.multAdd((1.0 / (boundFairness[1] - boundFairness[0])), fairnessSum);
        model.addConstr(fairness, GRB.EQUAL, fairnessFinal, "fairness");

        // Set solver parameters
        model.set(GRB.DoubleParam.TimeLimit, this.timeLimit); // Set time limit for the solver
        model.set(GRB.DoubleParam.MIPGap, 1e-3); // Set MIP gap tolerance for optimality

        //**********************************************************************
        //**********************************************************************

        // Start the optimization and measure computation time
        double startT = System.currentTimeMillis();
        model.optimize();
        double computationTime = (System.currentTimeMillis() - startT) / 1000.0;

        // Output the results after optimization
        System.out.println("########################## ff: " + model.get(GRB.DoubleAttr.ObjVal));
        System.out.println("########################## MIP GAP: \t" + model.get(GRB.DoubleAttr.MIPGap));
        System.out.println("########################## Fairness: \t" + fairness.get(GRB.DoubleAttr.X));

        // Arrays to store the results of optimization
        int[][][] resX = new int[instance.GreenType.length][instance.W][instance.H];
        int nCluster = 0;
        double[][][] resZ = new double[instance.UrbanChallenges.length][instance.W][instance.H];
        double[][][] resMu = new double[instance.UrbanChallenges.length][instance.W][instance.H];
        double[][][] resZbar = new double[instance.UrbanChallenges.length][instance.W][instance.H];
        double[][][] A_new = new double[instance.UrbanChallenges.length][instance.W][instance.H];

        // Extract the optimization results from the variables
        for (int i = 0; i < instance.W; i++) {
            for (int j = 0; j < instance.H; j++) {
                for (int b = 0; b < instance.UrbanChallenges.length; b++) {
                    // Store the results for Z, Mu, Z_bar, and updated A values
                    resZ[b][i][j] = z_var[b][i][j].get(GRB.DoubleAttr.X);
                    resMu[b][i][j] = z_mu[b][i][j].get(GRB.DoubleAttr.X);
                    resZbar[b][i][j] = z_bar_var[b][i][j].get(GRB.DoubleAttr.X);
                    A_new[b][i][j] = instance.A[b][i][j] - resZbar[b][i][j];
                }
                for (int t = 0; t < instance.GreenType.length; t++) {
                    // Store the result of whether a green type is assigned to the current tile
                    if (x_var[t][i][j].get(GRB.DoubleAttr.X) > 0.9) {
                        resX[t][i][j] = 1;
                    } else {
                        resX[t][i][j] = 0;
                    }
                }
            }
        }

        // Check if any clusters are used and count the number of active clusters
        for (int b = 0, l2 = this.clusterSize.size(); b < l2; b++) {
            if (clusterSize.get(b) > 0 && clusters[b].get(GRB.DoubleAttr.X) > 0) {
                nCluster++;
            }
        }

        // Denormalize the result data
        A_new = InstanceUtils.denormalizeData(instance, A_new);
        instance.A = InstanceUtils.denormalizeData(instance, instance.A);

        // Save the results (Objective value, MIP gap, computation time, etc.)
        ResultsUtils.saveResults(instance, config, model.get(GRB.DoubleAttr.ObjVal), model.get(GRB.DoubleAttr.MIPGap), computationTime, resX, resZ, A_new, kernelSet, mat_of_forbidden, nCluster);
    }
}
