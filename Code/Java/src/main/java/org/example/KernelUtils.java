package org.example;

import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

/**
 * Utility class for managing kernels.
 * Provides methods to read, print, and normalize kernels.
 */
public class KernelUtils {

    /***
     * Reads kernel data from a JSON file at the specified path and returns a nested HashMap
     * organizing kernels by "UrbanChallenge" and "GreenType".
     *
     * @param path Path to the JSON file containing kernel data.
     * @return A nested HashMap<String, HashMap<String, Kernel>> structure mapping UrbanChallenges to GreenTypes and their corresponding kernels.
     */
    public static HashMap<String, HashMap<String, Kernel>> readKernels(String path) {
        // Initialize the main map to store kernels.
        HashMap<String, HashMap<String, Kernel>> kernelSet = new HashMap<>();
        JSONParser jsonParser = new JSONParser();

        try (FileReader reader = new FileReader(path)) {
            // Parse the JSON file.
            JSONObject json = (JSONObject) jsonParser.parse(reader);

            // Extract the necessary arrays and objects from the JSON structure.
            JSONArray GreenType = (JSONArray) json.get(Config.GT);
            JSONArray UCs = (JSONArray) json.get(Config.UC);
            JSONObject SizeK = (JSONObject) json.get(Config.Sizes);
            JSONObject K = (JSONObject) json.get(Config.K);

            // Loop through each UrbanChallenge.
            for (int b = 0, l2 = UCs.size(); b < l2; b++) {
                String currB = UCs.get(b).toString();

                // Ensure the kernelSet contains an entry for the current UC.
                kernelSet.putIfAbsent(currB, new HashMap<>());

                // Retrieve size and kernel data for the current UC.
                JSONObject currSizeByB = (JSONObject) SizeK.get(currB);
                JSONObject currKByB = (JSONObject) K.get(currB);

                // Loop through each GreenType.
                for (int g = 0, l = GreenType.size(); g < l; g++) {
                    String currG = GreenType.get(g).toString();

                    // Extract kernel dimensions for the current GreenType.
                    int currW = Integer.parseInt(((JSONObject) currSizeByB.get(currG)).get(Config.row).toString());
                    int currH = Integer.parseInt(((JSONObject) currSizeByB.get(currG)).get(Config.col).toString());

                    // Read the kernel data into a 2D array.
                    JSONArray kToParse = (JSONArray) currKByB.get(currG);
                    double[][] currK = new double[currH][currW];
                    for (int w = 0; w < currW; w++) {
                        for (int h = 0; h < currH; h++) {
                            int refInVector = (w * currW) + h; // Calculate linear index.
                            currK[w][h] = Double.parseDouble(kToParse.get(refInVector).toString());
                        }
                    }

                    // Add the kernel to the nested map.
                    kernelSet.get(currB).put(currG, new Kernel(currG, currB, currW, currH, currK));
                }
            }
        } catch (IOException | ParseException e) {
            // Handle file I/O and parsing errors.
            e.printStackTrace();
        }

        return kernelSet; // Return the completed kernel set.
    }

    /***
     * Prints all kernels in the provided kernel set.
     *
     * @param kernelSet The kernel set to print.
     */
    public static void printKernels(HashMap<String, HashMap<String, Kernel>> kernelSet) {
        // Iterate over each UrbanChallenge and its associated GreenTypes.
        kernelSet.forEach((b, ByG) -> {
            ByG.forEach((g, k) -> {
                // Print each kernel.
                k.printKernel();
            });
        });
    }

    /***
     * Normalizes the values of kernels in the provided kernel set based on the given instance's
     * UrbanChallenge boundaries.
     *
     * @param instance The instance providing normalization boundaries.
     * @param kernelSet The kernel set to normalize.
     * @return A normalized kernel set.
     */
    public static HashMap<String, HashMap<String, Kernel>> normalizeKernel(Instance instance, HashMap<String, HashMap<String, Kernel>> kernelSet) {
        // Loop through each UC in the instance.
        for (int b = 0; b < instance.UrbanChallenges.length; b++) {
            String bS = instance.UrbanChallenges[b]; // Current UC as a string.

            // Loop through each GreenType in the instance.
            for (int g = 0; g < instance.GreenType.length; g++) {
                String gS = instance.GreenType[g]; // Current GreenType as a string.

                // Retrieve kernel dimensions.
                int wK = kernelSet.get(bS).get(gS).w;
                int hK = kernelSet.get(bS).get(gS).h;

                // Normalize each kernel value.
                for (int w = 0; w < wK; w++) {
                    for (int h = 0; h < hK; h++) {
                        // Apply normalization formula.
                        kernelSet.get(bS).get(gS).k[w][h] = (kernelSet.get(bS).get(gS).k[w][h] - instance.BoundUCs.get(bS)[0])
                                / (instance.BoundUCs.get(bS)[1] - instance.BoundUCs.get(bS)[0]);

                        // Correct out-of-range values.
                        if (kernelSet.get(bS).get(gS).k[w][h] > 1.0 || kernelSet.get(bS).get(gS).k[w][h] < 0.0 || Double.isNaN(kernelSet.get(bS).get(gS).k[w][h])) {
                            kernelSet.get(bS).get(gS).k[w][h] = kernelSet.get(bS).get(gS).k[w][h] > 1.0 ? 1.0 : 0.0;
                        }
                    }
                }
            }
        }
        return kernelSet; // Return the normalized kernel set.
    }
}
