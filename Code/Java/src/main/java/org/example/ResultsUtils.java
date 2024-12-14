package org.example;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Utility class for managing and saving results of an optimization problem.
 */
public class ResultsUtils {
    // Static fields to store general information about the instance and configuration
    static String name; // Name of the instance
    static String path; // Output path for saving results
    static int W; // Width of the grid
    static int H; // Height of the grid
    static HashMap<String, Double> initFF; // Initial objective function values

    /**
     * Initializes the utility with basic information about the instance and configuration.
     *
     * @param instance_name Name of the instance.
     * @param instance      The problem instance.
     * @param outPath       Output path for saving results.
     * @param singleFF_init Initial objective function values.
     */
    public static void getInfos(String instance_name, Instance instance, String outPath, HashMap<String, Double> singleFF_init) {
        name = instance_name;
        path = outPath;
        initFF = singleFF_init;
        W = instance.W;
        H = instance.H;
    }

    /**
     * Saves the results of the optimization process to a JSON file.
     *
     * @param instance         The problem instance.
     * @param config           The configuration of the optimization.
     * @param objValue         The final objective value.
     * @param mipGap           The final gap for mixed-integer programming.
     * @param time             Total time taken for the optimization.
     * @param resX             The solution matrix for variable X.
     * @param resZ             The solution matrix for variable Z.
     * @param A_new            Updated auxiliary matrix.
     * @param kernel           Kernel structure mapping.
     * @param mat_of_forbidden Matrix of forbidden zones.
     * @param nCluster         Number of clusters used in the optimization.
     */
    public static void saveResults(Instance instance, Config config, double objValue, double mipGap, double time, int[][][] resX, double[][][] resZ, double[][][] A_new, HashMap<String, HashMap<String, Kernel>> kernel, int[][] mat_of_forbidden, int nCluster) {
        JSONObject json = new JSONObject();

        // Add basic information about the instance and results
        json.put("W", instance.W);
        json.put("H", instance.H);
        json.put("GreenType", Arrays.stream(instance.GreenType).toList());
        json.put("UrbanChallenges", Arrays.stream(instance.UrbanChallenges).toList());
        json.put("Time", time);
        json.put("MipGap", mipGap);
        json.put("Of", objValue);
        json.put("Budget", config.Budget);
        json.put("#Cluster", nCluster);


        // Compute and add the final objective function values
        HashMap<String, Double> singleFF_final = InstanceUtils.computeOF(instance, config, A_new, mat_of_forbidden, kernel);
        JSONObject Of_list_final = new JSONObject();
        for (String key : singleFF_final.keySet()) {
            Of_list_final.put(key, singleFF_final.get(key));
        }
        json.put("Of_list_final", Of_list_final);

        // Add the initial objective function values
        JSONObject Of_list_init = new JSONObject();
        for (String key : initFF.keySet()) {
            Of_list_init.put(key, initFF.get(key));
        }
        json.put("Of_list_init", Of_list_init);

        // Compute fairness metrics using Gini coefficient
        double gini_init = InstanceUtils.computeGiniCoefficient(instance, generatePreExistentMatrix(instance), kernel.get(instance.UrbanChallenges[instance.UrbanChallenges.length - 1]));
        double gini_final = InstanceUtils.computeGiniCoefficient(instance, resX, kernel.get(instance.UrbanChallenges[instance.UrbanChallenges.length - 1]));

        if (gini_init < gini_final) {
            System.out.println("POSSIBLE ERROR IN GINI COMPUTATION: init=" + gini_init + " final=" + gini_final);
            gini_final = gini_init;
        }
        json.put("FairnessGini_init", gini_init);
        json.put("FairnessGini_final", gini_final);

        // Process results for installations and budget
        processInstallationsAndBudget(instance, config, resX, json);

        // Process results for coverage (variable Z)
        processCoverage(instance, A_new, json);

        // Write the JSON object to a file
        try (FileWriter file = new FileWriter(path + "Out_" + name)) {
            file.write(json.toString());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Generates a 3D matrix representing pre-existing installations.
     * The matrix has dimensions corresponding to the number of green types,
     * the width, and the height of the instance. Each element is marked with 1
     * if a pre-existing installation is present at the corresponding location.
     *
     * @param instance the Instance object containing the pre-existing installations
     * @return a 3D array representing the pre-existing installations matrix
     */
    private static int[][][] generatePreExistentMatrix(Instance instance) {
        int[][][] preExistentX = new int[instance.GreenType.length][instance.W][instance.H];

        // Iterate over each green type and set the corresponding cells with pre-existing installations
        for (int i = 0; i < instance.GreenType.length; i++) {
            for (int cellIndex = 0, size = instance.PreExistent.get(instance.GreenType[i]).size(); cellIndex < size; cellIndex++) {
                Cells_Tile cell = instance.PreExistent.get(instance.GreenType[i]).get(cellIndex);
                preExistentX[i][cell.r][cell.c] = 1; // Mark the cell with 1 if it has a pre-existing installation
            }
        }
        return preExistentX; // Return the matrix
    }

    /**
     * Processes the results related to installations (green types) and budget,
     * and updates the given JSON object with relevant information such as
     * the number of installations, coverage, and budget spent.
     *
     * @param instance the Instance object containing the green types and configuration
     * @param config   the configuration object containing cost information
     * @param resX     the 3D array representing the installation results (1 for installed)
     * @param json     the JSON object to be updated with the results
     */
    private static void processInstallationsAndBudget(Instance instance, Config config, int[][][] resX, JSONObject json) {
        double spentBudget = 0.0;
        int totalX = 0;

        // JSON objects to store the results
        JSONObject X = new JSONObject();
        JSONObject numberOfX = new JSONObject();
        JSONObject coveragePerNBS = new JSONObject();
        JSONObject BudgetSpentPerNBS = new JSONObject();

        // Process each green type
        for (int i = 0; i < instance.GreenType.length; i++) {
            int countCurrX = 0;
            JSONArray G_type = new JSONArray();

            // Iterate through each cell and record installations and counts
            for (int w = 0; w < resX[i].length; w++) {
                for (int h = 0; h < resX[i][w].length; h++) {
                    G_type.add(resX[i][w][h]); // Add the installation result (1 or 0) for the cell
                    if (resX[i][w][h] == 1 && !instance.checkPreExist(i, w, h)) {
                        countCurrX++; // Increment count for new installations (not pre-existing)
                    }
                }
            }

            // Update the JSON with data for this green type
            X.put(instance.GreenType[i], G_type);
            numberOfX.put(instance.GreenType[i], countCurrX);
            coveragePerNBS.put(instance.GreenType[i], (double) countCurrX / (instance.W * instance.H));
            BudgetSpentPerNBS.put(instance.GreenType[i], countCurrX * config.cost.get(instance.GreenType[i]));

            // Accumulate total installations and budget spent
            totalX += countCurrX;
            spentBudget += countCurrX * config.cost.get(instance.GreenType[i]);
        }

        // Update the JSON with the aggregated results
        json.put("X", X);
        json.put("InstalledNBS", numberOfX);
        json.put("TotalInstalledNBS", totalX);
        json.put("BudgetSpent", spentBudget);
        json.put("BudgetSpentPerNBS", BudgetSpentPerNBS);
        json.put("%CoverageTotal", (double) totalX / (instance.W * instance.H));
    }

    /**
     * Processes the coverage results (represented by variable Z) and updates
     * the JSON object with coverage data for each urban challenge.
     *
     * @param instance the Instance object containing the urban challenges
     * @param A_new    the 3D array of coverage data (Z) for each urban challenge
     * @param json     the JSON object to be updated with the coverage results
     */
    private static void processCoverage(Instance instance, double[][][] A_new, JSONObject json) {
        JSONObject Z = new JSONObject();

        // Process coverage for each urban challenge
        for (int b = 0; b < instance.UrbanChallenges.length; b++) {
            JSONArray B_type = new JSONArray();

            // Iterate through each cell and add the coverage value (Z) to the JSON
            for (int w = 0; w < A_new[b].length; w++) {
                for (int h = 0; h < A_new[b][w].length; h++) {
                    B_type.add(A_new[b][w][h]); // Add the coverage value for this cell
                }
            }

            // Update the JSON with the coverage data for the current urban challenge
            Z.put(instance.UrbanChallenges[b], B_type);
        }

        // Add the coverage data to the JSON
        json.put("Z", Z);
    }
}
