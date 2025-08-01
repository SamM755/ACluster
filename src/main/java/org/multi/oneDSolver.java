package org.multi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class oneDSolver {

    
    private static class IndexedValue implements Comparable<IndexedValue> {
        long value;
        int originalIndex;

        IndexedValue(long value, int originalIndex) {
            this.value = value;
            this.originalIndex = originalIndex;
        }

        @Override
        public int compareTo(IndexedValue other) {
            return Long.compare(this.value, other.value);
        }
    }

    

    public static Object[] optimalKMedoid1D(long[][] data, int k, int dim) {
        int n = data.length;
        if (n == 0 || k <= 0) {
            return new Object[]{new long[0][1], new int[0], new int[0]};
        }
        if (k >= n) {
            
            long[][] medoids = new long[n][1];
            int[] assignments = new int[n];
            int[] sizes = new int[n];
            for (int i = 0; i < n; i++) {
                medoids[i] = new long[]{data[i][0]};
                assignments[i] = i;
                sizes[i] = 1;
            }
            return sortResults(medoids, assignments, sizes);
        }

        
        IndexedValue[] sortedData = new IndexedValue[n];
        for (int i = 0; i < n; i++) {
            sortedData[i] = new IndexedValue(data[i][0], i);
        }
        Arrays.sort(sortedData);

        
        long[][] costs = precomputeCosts(sortedData);

        
        
        long[][] dp = new long[k + 1][n + 1];
        
        int[][] parent = new int[k + 1][n + 1];

        
        for (int i = 0; i <= k; i++) {
            Arrays.fill(dp[i], Long.MAX_VALUE / 2);
        }
        dp[0][0] = 0; 

        
        for (int i = 1; i <= k; i++) {
            
            knuthYaoCompute(i, 1, n, 0, n - 1, dp, parent, costs);
        }

        
        List<long[]> medoidsList = new ArrayList<>();
        int[] assignments = new int[n];
        List<Integer> clusterSizes = new ArrayList<>();

        int currentN = n; 
        for (int currentK = k; currentK > 0; currentK--) {
            int p = parent[currentK][currentN]; 
            int start_idx_sorted = p;
            int end_idx_sorted = currentN - 1;

            if (p < 0 || start_idx_sorted > end_idx_sorted) continue; 

            
            int medianIndexInSorted = start_idx_sorted + (end_idx_sorted - start_idx_sorted) / 2;
            medoidsList.add(new long[]{sortedData[medianIndexInSorted].value});

            
            for (int i = start_idx_sorted; i <= end_idx_sorted; i++) {
                assignments[sortedData[i].originalIndex] = currentK - 1; 
            }
            clusterSizes.add(end_idx_sorted - start_idx_sorted + 1);

            
            currentN = p;
        }

        
        Collections.reverse(medoidsList);
        Collections.reverse(clusterSizes);

        
        long[][] finalMedoids = medoidsList.toArray(new long[0][]);
        int[] finalSizes = clusterSizes.stream().mapToInt(i -> i).toArray();

        
        return sortResults(finalMedoids, assignments, finalSizes);
    }

    private static void knuthYaoCompute(int i, int j_low, int j_high, int p_low, int p_high,
                                        long[][] dp, int[][] parent, long[][] costs) {
        if (j_low > j_high) {
            return;
        }

        int j_mid = j_low + (j_high - j_low) / 2;
        long min_cost = Long.MAX_VALUE;

        
        
        
        
        
        int p_best = p_low;
        

        
        
        int p_end = Math.min(j_mid - 1, p_high);

        for (int p = p_low; p <= p_end; p++) {
            
            
            
            if (dp[i-1][p] < Long.MAX_VALUE / 2) {
                long current_cost = dp[i-1][p] + costs[p][j_mid - 1];
                if (current_cost < min_cost) {
                    min_cost = current_cost;
                    p_best = p; 
                }
            }
        }

        dp[i][j_mid] = min_cost;
        parent[i][j_mid] = p_best;

        
        
        knuthYaoCompute(i, j_low, j_mid - 1, p_low, p_best, dp, parent, costs);
        
        knuthYaoCompute(i, j_mid + 1, j_high, p_best, p_high, dp, parent, costs);
    }

    public static Object[] optimalAdaptiveBPRE1D(long[][] pageData, int dim) {
        int n = pageData.length;
        dim = 1; 
        if (n == 0) {
            return new Object[]{new long[0][1], new int[0], new int[0]};
        }

        final long BASE_POINT_COST = (long) dim * 64;

        
        IndexedValue[] sortedData = new IndexedValue[n];
        for (int i = 0; i < n; i++) {
            sortedData[i] = new IndexedValue(pageData[i][0], i);
        }
        Arrays.sort(sortedData);

        
        long[][] costs = precomputeCosts(sortedData);

        
        
        
        long[] dp = new long[n + 1];
        int[] parent = new int[n + 1];

        Arrays.fill(dp, Long.MAX_VALUE / 2);
        dp[0] = 0;

        
        
        compute(1, n, 0, n - 1, dp, parent, costs, BASE_POINT_COST);
        

        
        List<long[]> medoidsList = new ArrayList<>();
        int[] assignments = new int[n];
        List<Integer> clusterSizes = new ArrayList<>();

        int currentN = n;
        int clusterId = 0;
        while (currentN > 0) {
            int start = parent[currentN];
            int end = currentN - 1;

            int medianIndexInSorted = start + (end - start) / 2;
            medoidsList.add(new long[]{sortedData[medianIndexInSorted].value});

            for (int i = start; i <= end; i++) {
                assignments[sortedData[i].originalIndex] = clusterId;
            }
            clusterSizes.add(end - start + 1);

            currentN = start;
            clusterId++;
        }

        Collections.reverse(medoidsList);
        Collections.reverse(clusterSizes);

        int numClusters = medoidsList.size();
        int[] remappedAssignments = new int[n];
        for (int i = 0; i < n; i++) {
            remappedAssignments[i] = numClusters - 1 - assignments[i];
        }

        
        long[][] finalMedoids = medoidsList.toArray(new long[0][]);
        int[] finalSizes = clusterSizes.stream().mapToInt(i -> i).toArray();

        return sortResults(finalMedoids, remappedAssignments, finalSizes);
    }

    /**
     * Divide and Conquer DP optimization helper function.
     * @param i_low         The lower bound of the dp array indices to compute.
     * @param i_high        The upper bound of the dp array indices to compute.
     * @param j_low         The lower bound for the optimal split point 'j'.
     * @param j_high        The upper bound for the optimal split point 'j'.
     * @param dp            The DP table.
     * @param parent        The parent pointer table for backtracking.
     * @param costs         The precomputed cost matrix.
     * @param basePointCost The fixed cost for adding a new base point.
     */
    private static void compute(int i_low, int i_high, int j_low, int j_high,
                                long[] dp, int[] parent, long[][] costs, final long basePointCost) {
        if (i_low > i_high) {
            return;
        }

        int i_mid = i_low + (i_high - i_low) / 2;
        long min_cost = Long.MAX_VALUE;
        int j_opt = -1;

        
        
        int j_end = Math.min(i_mid - 1, j_high);
        for (int j = j_low; j <= j_end; j++) {
            
            
            if (dp[j] < Long.MAX_VALUE / 2) {
                long currentCost = dp[j] + costs[j][i_mid - 1] + basePointCost;
                if (currentCost < min_cost) {
                    min_cost = currentCost;
                    j_opt = j;
                }
            }
        }

        dp[i_mid] = min_cost;
        parent[i_mid] = j_opt;

        
        
        compute(i_low, i_mid - 1, j_low, j_opt, dp, parent, costs, basePointCost);
        
        compute(i_mid + 1, i_high, j_opt, j_high, dp, parent, costs, basePointCost);
    }

    
    private static long[][] precomputeCosts(IndexedValue[] sortedData) {
        int n = sortedData.length;
        long[][] costs = new long[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                
                int medianIndex = i + (j - i) / 2;
                long medianValue = sortedData[medianIndex].value;
                
                long sumDist = 0;
                for (int l = i; l <= j; l++) {
                    sumDist += Math.abs(sortedData[l].value - medianValue);
                }
                costs[i][j] = sumDist;
            }
        }
        return costs;
    }

    
    private static Object[] sortResults(long[][] medoids, int[] clusterAssignment, int[] clusterSizes) {
        int k = medoids.length;
        if (k == 0) {
            return new Object[]{new long[0][1], new int[0], new int[0]};
        }

        
        class ClusterInfo implements Comparable<ClusterInfo> {
            long[] medoid;
            int size;
            int originalIndex;

            ClusterInfo(long[] medoid, int size, int originalIndex) {
                this.medoid = medoid;
                this.size = size;
                this.originalIndex = originalIndex;
            }

            @Override
            public int compareTo(ClusterInfo other) {
                return Integer.compare(other.size, this.size); 
            }
        }

        List<ClusterInfo> clusterInfoList = new ArrayList<>();
        for (int i = 0; i < k; i++) {
            clusterInfoList.add(new ClusterInfo(medoids[i], clusterSizes[i], i));
        }

        
        Collections.sort(clusterInfoList);

        
        int[] indexMap = new int[k];
        for (int i = 0; i < k; i++) {
            indexMap[clusterInfoList.get(i).originalIndex] = i;
        }

        
        long[][] sortedMedoids = new long[k][medoids[0].length];
        long[] sortedClusterSizes = new long[k];
        for (int i = 0; i < k; i++) {
            sortedMedoids[i] = clusterInfoList.get(i).medoid;
            sortedClusterSizes[i] = clusterInfoList.get(i).size;
        }

        int[] newClusterAssignment = new int[clusterAssignment.length];
        for (int i = 0; i < clusterAssignment.length; i++) {
            newClusterAssignment[i] = indexMap[clusterAssignment[i]];
        }

        return new Object[]{sortedMedoids, newClusterAssignment, sortedClusterSizes};
    }

    
    public static void main(String[] args) {
        
        long[][] data1D = {
                {10}, {12}, {15}, {30}, {32}, {33}, {50}, {55}, {58}, {100}
        };
        int n = data1D.length;

        int k = 3;
        Object[] resultKMedoid = optimalKMedoid1D(data1D, k, 1);

        long[][] medoidsK = (long[][]) resultKMedoid[0];
        int[] assignmentsK = (int[]) resultKMedoid[1];
        int[] sizesK = (int[]) resultKMedoid[2];

        for (long[] medoid : medoidsK) {
            System.out.println(Arrays.toString(medoid));
        }

        Object[] resultABPRE = optimalAdaptiveBPRE1D(data1D, 1);
        long[][] medoidsA = (long[][]) resultABPRE[0];
        int[] assignmentsA = (int[]) resultABPRE[1];
        int[] sizesA = (int[]) resultABPRE[2];

        for (long[] medoid : medoidsA) {
            System.out.println(Arrays.toString(medoid));
        }
        
    }
}
