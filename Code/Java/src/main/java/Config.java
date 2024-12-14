import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

public class Config {
    // FOLDER CONFIGURATION: Paths for various folders used in the program
    public static String KernelFolder = "Config/"; // Folder containing kernel data
    public static String KernelName = "Kernels.json"; // Name of the kernel file
    public static String TestInstanceFolder = "Test/"; // Folder containing test instances
    public static String OutputFolder = "Results/"; // Folder where results will be saved

    // MILP (Mixed Integer Linear Programming) PARAMETERS
    public static int timeLimit = 1800; // Time limit for the MILP solver (in seconds)
    public static boolean showOut = false; // Whether to show output or not
    public static double[] theta = { // Array of theta values used in the MILP model
            0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    };

    // CLUSTER CREATION PARAMETERS: Define properties for clustering
    public static int maxClusterSize = 50; // Maximum size for each cluster
    public static int maxLimitSize = (int) (Config.maxClusterSize * 1.5); // Maximum limit for cluster size (1.5x maxClusterSize)
    public static int minClusterSize = 5; // Minimum size for each cluster
    public static String ClusteredGT = "UrbanPark"; // Type of green space for clustering (Urban Park)

    // BUDGET PARAMETERS: Budget range parameters for the solution
    public static double percBudgetMax = 0.5; // Maximum percentage of budget (0.5 means 50% of max cost)
    public static double percBudgetMin = 0.3; // Minimum percentage of budget (0.3 means 30% of max cost)

    // URBAN CHALLENGE (UC) PARAMETERS: Parameter for reducing urban challenges
    public static double percBoundUCs = 0.2; // Percentage reduction allowed for urban challenges (1 means no reduction, 0 means full reduction)

    // STRUCTURE FOR HANDLING INSTANCES
    public HashMap<String, Double> cost = new HashMap<>(); // Maps the green type to its cost
    public HashMap<String, Double> normalizedCost = new HashMap<>(); // Normalized cost map
    public HashMap<String, String[]> UCMap = new HashMap<>(); // Maps urban challenges to their objectives (Avg, Max, Sum, Gini)
    public ArrayList<String> UCsOF = new ArrayList<>(); // List of objective functions for urban challenges
    public double Budget; // Budget for the instance
    public double NormalizedBudget; // Normalized budget value
    public double[] BigM; // Array for BigM values, used for defining bounds in MILP models

    // KERNEL STRINGS: Keys used for identifying kernel-related data
    public static String K = "K"; // Kernel matrix
    public static String Sizes = "SizeK"; // Sizes for kernels
    public static String row = "r"; // Row identifier in kernel matrix
    public static String col = "c"; // Column identifier in kernel matrix

    // INSTANCE STRINGS: Keys for instance-related data
    public static String W = "W"; // Width of the instance
    public static String H = "H"; // Height of the instance
    public static String GT = "GreenType"; // Type of green space
    public static String UC = "UrbanChallenge"; // Urban challenges to address
    public static String PreExist = "PreExistent"; // Pre-existing elements in the instance
    public static String Forbid = "Forbidden"; // Forbidden elements in the instance
    public static String A = "A"; // A matrix related to the instance data
    public static String[] dataName = { // Different types of data available for the instance
            "TempMax", "TempMin", "Pm10", "Pm2", "Fairness"
    };

    // MODEL STRINGS: Defines different objectives in the MILP model
    public static String Avg = "Avg"; // Average value
    public static String Max = "Max"; // Maximum value

    // Constructor for the Config class
    public Config() {
        // Initialize NBS (Nature-based Solutions) costs
        this.cost = new HashMap<>();
        this.cost.put("GreenWall", 780.9); // Cost for Green Wall
        this.cost.put("GreenRoof", 520.0); // Cost for Green Roof
        this.cost.put("StreetTree", 210.0); // Cost for Street Tree
        this.cost.put("UrbanPark", 370.8); // Cost for Urban Park

        this.normalizedCost = new HashMap<>(); // Initialize the normalized cost map

        // Define the objective functions (OF) for different urban challenges
        this.UCMap.put("TempMax", new String[]{"Avg", "Max"});
        this.UCMap.put("TempMin", new String[]{"Avg", "Max"});
        this.UCMap.put("Pm10", new String[]{"Avg", "Max"});
        this.UCMap.put("Pm2", new String[]{"Avg", "Max"});
        this.UCMap.put("Costs", new String[]{"Sum"});
        this.UCMap.put("Fairness", new String[]{"Gini"});

        // Add the combination of challenges and objectives to UCsOF
        this.UCsOF = new ArrayList<>();
        for (String k : this.UCMap.keySet()) {
            for (String s : this.UCMap.get(k)) {
                UCsOF.add(k + s);
            }
        }

        // Set the BigM values for different data types
        this.BigM = new double[this.dataName.length];
        this.BigM[0] = 10.0; // Max value for TempMax
        this.BigM[1] = 10.0; // Max value for TempMin
        this.BigM[2] = 10.0; // Max value for Pm10
        this.BigM[3] = 10.0; // Max value for Pm2
        this.BigM[4] = 100000.0; // Max value for Fairness (large value for fairness-related constraints)
    }

    // Method to compute the maximum budget based on the width (W) and height (H) of the instance
    public void budgetComputationMax(int W, int H) {
        double max = 0;
        // Find the maximum cost for the different green types
        for (String key : cost.keySet()) {
            if (cost.get(key) > max) {
                max = cost.get(key); // Update the max if the current cost is higher
            }
        }

        // Generate a random budget percentage between the minimum and maximum budget percentage
        Random random = new Random();
        double percBudget = random.nextDouble(percBudgetMax - percBudgetMin) + percBudgetMin;
        // Compute the budget as a product of the width, height, max cost, and the random percentage
        this.Budget = W * H * max * percBudget;
    }
}
