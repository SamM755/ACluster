package org.multi;

import java.lang.reflect.Member;
import java.util.*;
import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.VectorSpecies;


class MedoidSortHelper implements Comparable<MedoidSortHelper> {
    long[] medoid;
    long size;
    int originalIndex;


    MedoidSortHelper(long[] medoid, long size, int originalIndex) {
        this.medoid = medoid;
        this.size = size;
        this.originalIndex = originalIndex;
    }

    @Override
    public int compareTo(MedoidSortHelper other) {
        return Long.compare(this.size, other.size);
    }
}

public class Base {
    private static final Random rand = new Random();
    private static final VectorSpecies<Long> SPECIES = LongVector.SPECIES_PREFERRED;
    private static final int LANE_COUNT = SPECIES.length();

    public static long manhattanDistance(long[] p1, long[] p2) {
        long sum = 0;
        for (int i = 0; i < p1.length; i++) {
            sum += Math.abs(p1[i] - p2[i]);
        }
        return sum;
    }

    public static long logEncodingCost(long value) {
        return (long) (Math.log(value + 1) / Math.log(2));
    }













    private static long bitLength(long value) {
        if(value==0) return 1;
        return 64 - Long.numberOfLeadingZeros(value);
    }









    public static long calculateTotalLogCost(long[] p1, long[] p2, int dim) {
        if (p1 == null || p2 == null || p1.length != dim || p2.length != dim) {
            throw new IllegalArgumentException("Invalid data points or dimension.");
        }
        long totalCost = 0;
        for (int i = 0; i < dim; i++) {
            long residual = Math.abs(p1[i] - p2[i]);
            totalCost += logEncodingCost(residual) + 1;

        }
        return totalCost;
    }

    private static long calculateBasePointStorageCost_simd(long[] point, int dim) {
        long totalBits = 0;
        for (int i = 0; i < dim; i++) {
            totalBits += 1 + bitLength(Math.abs(point[i]));
        }
        return totalBits;
    }

    public static long calculateTotalLogCost_Ultimate(long[] p1, long[] p2, int dim) {
        long totalCost = 0;
        int i = 0;

        for (; i < SPECIES.loopBound(dim); i += LANE_COUNT) {
            LongVector v1 = LongVector.fromArray(SPECIES, p1, i);
            LongVector v2 = LongVector.fromArray(SPECIES, p2, i);
            LongVector vResidual = v1.sub(v2).abs();

            for (int j = 0; j < LANE_COUNT; j++) {
                long residual = vResidual.lane(j);
                totalCost += 1 + bitLength(residual);
            }
        }

        for (; i < dim; i++) {
            long residual = Math.abs(p1[i] - p2[i]);
            totalCost += 1 + bitLength(residual);
        }

        return totalCost;
    }

    public static long[][] updateMedoids(long[][] data, int[] clusters, long[][] currentMedoids, int dim) {
        long[][] newMedoids = new long[currentMedoids.length][dim];
        Map<Integer, List<long[]>> clusterMap = new HashMap<>();

        for (int i = 0; i < data.length; i++) {
            clusterMap.computeIfAbsent(clusters[i], k -> new ArrayList<>()).add(data[i]);
        }

        for (int medoidIndex = 0; medoidIndex < currentMedoids.length; medoidIndex++) {
            List<long[]> clusterPoints = clusterMap.getOrDefault(medoidIndex, Collections.emptyList());
            if (clusterPoints.isEmpty()) {
                newMedoids[medoidIndex] = currentMedoids[medoidIndex];
                continue;
            }

            int m = clusterPoints.size();
            long[] totalCosts = new long[m];
            Arrays.fill(totalCosts, 0);

            for (int i = 0; i < m; i++) {
                long[] pointI = clusterPoints.get(i);
                for (int j = i + 1; j < m; j++) {
                    long[] pointJ = clusterPoints.get(j);
                    long pairCost = 0;
                    for (int d = 0; d < dim; d++) {
                        pairCost += logEncodingCost(Math.abs(pointI[d] - pointJ[d]));
                    }

                    

                    totalCosts[i] += pairCost;
                    totalCosts[j] += pairCost;
                }
            }

            int bestIndex = 0;
            long bestCost = totalCosts[0];
            for (int i = 1; i < m; i++) {
                if (totalCosts[i] < bestCost) {
                    bestCost = totalCosts[i];
                    bestIndex = i;
                }
            }
            newMedoids[medoidIndex] = clusterPoints.get(bestIndex);
        }
        return newMedoids;
    }


    private static int binarySearch(long[] prefixSums, long target) {
        int low = 0, high = prefixSums.length - 1;
        while (low < high) {
            int mid = (low + high) >>> 1;
            if (prefixSums[mid] >= target) {
                high = mid;
            } else {
                low = mid + 1;
            }
        }
        return low;
    }

    public static long[][] acceleratedInitialization(long[][] data, int k, int dim) {
        long[][] medoids = new long[k][dim];
        Set<String> selectedMedoids = new HashSet<>();

        int firstIndex = rand.nextInt(data.length);
        medoids[0] = data[firstIndex];
        selectedMedoids.add(Arrays.toString(medoids[0]));

        long[] distances = new long[data.length];
        for (int i = 0; i < data.length; i++) {
            distances[i] = manhattanDistance(data[i], data[firstIndex]);
        }

        




        for (int i = 1; i < k; i++) {

            long[] prefixSums = new long[data.length];
            prefixSums[0] = distances[0];
            for (int p = 1; p < data.length; p++) {
                prefixSums[p] = prefixSums[p - 1] + distances[p];
            }
            long totalDistance = prefixSums[data.length - 1];
            if (totalDistance == 0) {
                int idx = rand.nextInt(data.length);
                medoids[i] = data[idx];
                selectedMedoids.add(Arrays.toString(medoids[i]));
                continue;
            }

            long randValue = (long) (rand.nextDouble() * totalDistance);

            int chosenIdx = binarySearch(prefixSums, randValue);
            while (selectedMedoids.contains(Arrays.toString(data[chosenIdx]))) {
                chosenIdx = (chosenIdx + 1) % data.length;
            }

            medoids[i] = data[chosenIdx];
            selectedMedoids.add(Arrays.toString(medoids[i]));

            for (int idx = 0; idx < data.length; idx++) {
                long distNewMedoid = manhattanDistance(data[idx], medoids[i]);
                if (distNewMedoid < distances[idx]) {
                    distances[idx] = distNewMedoid;
                }
            }

            







        }

        return medoids;
    }

    public static Object[] kMedoidLogCost_SIMD(long[][] data, int k, int maxIter, double tol, int dim) {
        int n = data.length;

        Set<String> uniqueRows = new HashSet<>();
        for (int i = 0; i < data.length; i++) {
            uniqueRows.add(Arrays.toString(data[i]));
            if(uniqueRows.size()>k){
                break;
            }
        }
        int distinctCount = uniqueRows.size();
        if (distinctCount < k) {
            System.out.println("Distinct data points (" + distinctCount +
                    ") is less than input k (" + k + "), setting k to " + distinctCount);
            k = distinctCount;
        }

        long[][] medoids = acceleratedInitialization(data, k, dim);
        int[] clusterAssignment = new int[n];

        long previousTotalCost = Long.MAX_VALUE;
        for (int iteration = 0; iteration < maxIter; iteration++) {

            long totalCostThisRound = 0L;
            for (int i = 0; i < n; i++) {
                long minCost = Long.MAX_VALUE;
                int assignedMedoidIndex = -1;

                for (int m = 0; m < k; m++) {
                    long cost = calculateTotalLogCost(data[i], medoids[m], dim);

                    if (cost < minCost) {
                        minCost = cost;
                        assignedMedoidIndex = m;
                        if (minCost == 0) {
                            break;
                        }
                    }
                }
                clusterAssignment[i] = assignedMedoidIndex;
                totalCostThisRound += minCost;
            }

            if (iteration > 0 && (previousTotalCost - totalCostThisRound) < tol) {
                break;
            }
            previousTotalCost = totalCostThisRound;

            long[][] newMedoids = updateMedoids_SIMD(data, clusterAssignment, k, dim);

            if (Arrays.deepEquals(medoids, newMedoids)) {
                break;
            }
            medoids = newMedoids;
        }

        long[] finalClusterSizes = new long[k];
        for (int assignment : clusterAssignment) {
            if (assignment != -1) finalClusterSizes[assignment]++;
        }
        return sortResults(medoids, clusterAssignment, finalClusterSizes);
    }


    private static long[][] updateMedoids_SIMD(long[][] data, int[] clusterAssignment, int k, int dim) {
        long[][] newMedoids = new long[k][dim];
        List<Integer>[] clusterPoints = new ArrayList[k];
        for (int i = 0; i < k; i++) {
            clusterPoints[i] = new ArrayList<>();
        }
        for (int i = 0; i < data.length; i++) {
            if (clusterAssignment[i] != -1) {
                clusterPoints[clusterAssignment[i]].add(i);
            }
        }

        for (int m = 0; m < k; m++) {
            List<Integer> members = clusterPoints[m];
            if (members.isEmpty()) {
                continue;
            }

            long minTotalClusterCost = Long.MAX_VALUE;
            int newMedoidIndex = -1;

            for (Integer candidateIndex : members) {
                long currentCandidateTotalCost = 0L;
                for (Integer otherMemberIndex : members) {
                    currentCandidateTotalCost += calculateTotalLogCost(data[candidateIndex], data[otherMemberIndex], dim);
                }

                if (currentCandidateTotalCost < minTotalClusterCost) {
                    minTotalClusterCost = currentCandidateTotalCost;
                    newMedoidIndex = candidateIndex;
                }
            }
            if (newMedoidIndex != -1) {
                newMedoids[m] = data[newMedoidIndex];
            }
        }
        return newMedoids;
    }

    public static Object[] kMedoidLogCost(long[][] data, int k, int maxIter, double tol, int dim) {


        Set<String> uniqueRows = new HashSet<>();
        for (int i = 0; i < data.length; i++) {
            uniqueRows.add(Arrays.toString(data[i]));
            if(uniqueRows.size()>k){
                break;
            }
        }
        int distinctCount = uniqueRows.size();
        if (distinctCount < k) {
            System.out.println("Distinct data points (" + distinctCount +
                    ") is less than input k (" + k + "), setting k to " + distinctCount);
            k = distinctCount;
        }

        long[][] medoids = acceleratedInitialization(data, k, dim);


        int[] clusterAssignment = new int[data.length];
        long[] clusterSizes = new long[k];

        long previousTotalCost = Long.MAX_VALUE;
        for (int iteration = 0; iteration < maxIter; iteration++) {
            Arrays.fill(clusterSizes, 0);

            long totalCostThisRound = 0L;

            for (int i = 0; i < data.length; i++) {
                long minCost = Long.MAX_VALUE;
                int assignedMedoidIndex = -1;

                for (int m = 0; m < medoids.length; m++) {
                    boolean allZero = true;
                    long cost = 0;

                    for (int d = 0; d < dim; d++) {
                        long diff = data[i][d] - medoids[m][d];
                        if (diff != 0) {
                            allZero = false;
                            long absDiff = (diff > 0) ? diff : -diff;
                            cost += fastLogEncodingCost(absDiff);
                            if (cost > minCost) {
                                break;
                            }
                        }
                    }
                    if (allZero) {
                        minCost = 0;
                        assignedMedoidIndex = m;
                        break;
                    } else if (cost < minCost) {
                        minCost = cost;
                        assignedMedoidIndex = m;
                    }
                }

                clusterAssignment[i] = assignedMedoidIndex;
                clusterSizes[assignedMedoidIndex]++;
                totalCostThisRound += minCost;
            }

            double averageCost = (double) totalCostThisRound / data.length;



            if (iteration > 0 && (previousTotalCost - totalCostThisRound) < tol) {

                break;
            }
            previousTotalCost = totalCostThisRound;

            long[][] newMedoids = updateMedoids(data, clusterAssignment, medoids, dim);

            if (Arrays.deepEquals(medoids, newMedoids)) {

                break;
            }
            medoids = newMedoids;
        }

        Integer[] indices = new Integer[k];
        for (int i = 0; i < k; i++) indices[i] = i;
        Arrays.sort(indices, Comparator.comparingLong(i -> clusterSizes[i]));

        long[][] sortedMedoids = new long[k][dim];
        long[] sortedClusterSizes = new long[k];
        for (int i = 0; i < k; i++) {
            sortedMedoids[i] = medoids[indices[i]];
            sortedClusterSizes[i] = clusterSizes[indices[i]];
        }

        int[] newClusterAssignment = new int[data.length];
        List<Integer> indexList = Arrays.asList(indices);
        for (int i = 0; i < data.length; i++) {
            newClusterAssignment[i] = indexList.indexOf(clusterAssignment[i]);
        }

        return new Object[]{sortedMedoids, newClusterAssignment, sortedClusterSizes};
    }

    static class DataPointWithCost implements Comparable<DataPointWithCost> {
        int pointIndex;
        long residualCost;

        public DataPointWithCost(int pointIndex, long residualCost) {
            this.pointIndex = pointIndex;
            this.residualCost = residualCost;
        }

        
        
        @Override
        public int compareTo(DataPointWithCost other) {
            return Long.compare(this.residualCost, other.residualCost);
        }

        @Override
        public String toString() {
            return "Point(idx=" + pointIndex + ", cost=" + residualCost + ")";
        }
    }

    static class MedoidSortHelperOp implements Comparable<MedoidSortHelperOp> {
        long[] medoid;
        int originalIndex;
        long clusterSize;
        
        List<DataPointWithCost> points;

        public MedoidSortHelperOp(long[] medoid, int originalIndex, List<DataPointWithCost> points) {
            this.medoid = medoid;
            this.originalIndex = originalIndex;
            this.points = points;
            this.clusterSize = points.size();
        }

        @Override
        public int compareTo(MedoidSortHelperOp other) {
            
            return Long.compare(other.clusterSize, this.clusterSize);
        }
    }


    private static Object[] sortResults(long[][] medoids, int[] clusterAssignment, long[] clusterSize) {
        int k = medoids.length;
        List<MedoidSortHelper> sorters = new ArrayList<>();
        for (int i = 0; i < k; i++) {
            sorters.add(new MedoidSortHelper(medoids[i], clusterSize[i], i));
        }

        Collections.sort(sorters);

        long[][] sortedMedoids = new long[k][];
        long[] sortedClusterSize = new long[k];
        int[] oldToNewIndexMap = new int[k];

        for (int i = 0; i < k; i++) {
            MedoidSortHelper sortedItem = sorters.get(i);
            sortedMedoids[i] = sortedItem.medoid;
            sortedClusterSize[i] = sortedItem.size;
            oldToNewIndexMap[sortedItem.originalIndex] = i;
        }

        int[] sortedClusterAssignment = new int[clusterAssignment.length];
        for (int i = 0; i < clusterAssignment.length; i++) {
            int oldIndex = clusterAssignment[i];
            sortedClusterAssignment[i] = oldToNewIndexMap[oldIndex];
        }

        return new Object[]{sortedMedoids, sortedClusterAssignment, sortedClusterSize};
    }


    private static long fastLogEncodingCost(long diff) {

        if (diff < 0) {
            throw new IllegalArgumentException("diff < 0? unexpected: " + diff);
        }

        if (diff == 0) {
            return 0;
        }

        double val = Math.log1p(diff);
        return (long) Math.floor(val);

    }

    public static long calculateBasePointStorageCost(long[] basePoint, int dim) {
        if (basePoint == null || basePoint.length != dim) {
            throw new IllegalArgumentException("Invalid base point or dimension.");
        }
        long storageCost = 0;
        for (int i = 0; i < dim; i++) {
            storageCost += logEncodingCost(Math.abs(basePoint[i])) + 1;
        }
        return storageCost;
    }

    private static class ArrayWrapper {
        private final long[] data;
        private final int hashCode;

        public ArrayWrapper(long[] data) {
            this.data = data;
            this.hashCode = Arrays.hashCode(data);
        }

        @Override
        public int hashCode() {
            return this.hashCode;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            ArrayWrapper that = (ArrayWrapper) o;
            return Arrays.equals(this.data, that.data);
        }
    }

    public static Object[] adaptiveGreedyBasePointSelection(long[][] pageData, int dim) {



        int n = pageData.length;
        if (n == 0) {
            return new Object[]{new long[0][dim], new int[0], new long[0]};
        }


        List<long[]> medoidsList = new ArrayList<>();
        Set<ArrayWrapper> existingMedoidsSet = new HashSet<>();
        List<Set<Integer>> pointsInClusterList = new ArrayList<>();
        int[] pointToLeaderMap = new int[n];

        long[] firstPoint = Arrays.copyOf(pageData[0], dim);
        medoidsList.add(firstPoint);
        existingMedoidsSet.add(new ArrayWrapper(firstPoint));

        Set<Integer> firstClusterPoints = new HashSet<>();
        firstClusterPoints.add(0);
        pointsInClusterList.add(firstClusterPoints);
        pointToLeaderMap[0] = 0;

        for (int i = 1; i < n; i++) {
            long[] currentPoint = pageData[i];
            if (existingMedoidsSet.contains(new ArrayWrapper(currentPoint))) {
                int existingLeaderIndex = -1;
                for (int j = 0; j < medoidsList.size(); j++) {
                    if (Arrays.equals(medoidsList.get(j), currentPoint)) {
                        existingLeaderIndex = j;
                        break;
                    }
                }
                if (existingLeaderIndex != -1) {
                    pointsInClusterList.get(existingLeaderIndex).add(i);
                    pointToLeaderMap[i] = existingLeaderIndex;
                }
                continue;
            }

            
            int bestLeaderIndex = -1;
            long minCostToExistingLeader = Long.MAX_VALUE;
            for (int j = 0; j < medoidsList.size(); j++) {
                long cost = calculateTotalLogCost(currentPoint, medoidsList.get(j), dim);
                if (cost < minCostToExistingLeader) {
                    minCostToExistingLeader = cost;
                    bestLeaderIndex = j;
                }
            }

            long savingsFromCurrentPoint = minCostToExistingLeader;
            long savingsFromReassignment = 0;

            Set<Integer> pointsInBestLeaderCluster = pointsInClusterList.get(bestLeaderIndex);
            for (int pointIndexInCluster : pointsInBestLeaderCluster) {
                long[] p = pageData[pointIndexInCluster];
                long costToOldLeader = calculateTotalLogCost(p, medoidsList.get(bestLeaderIndex), dim);
                long costToNewPotentialLeader = calculateTotalLogCost(p, currentPoint, dim);
                if (costToNewPotentialLeader < costToOldLeader) {
                    savingsFromReassignment += (costToOldLeader - costToNewPotentialLeader);
                }
            }
            long totalSavings = savingsFromCurrentPoint + savingsFromReassignment;

            long storageCostForNewPoint = calculateBasePointStorageCost(currentPoint, dim);
            if (totalSavings > storageCostForNewPoint) {
                int newLeaderId = medoidsList.size();
                medoidsList.add(Arrays.copyOf(currentPoint, dim));
                existingMedoidsSet.add(new ArrayWrapper(currentPoint));

                Set<Integer> newClusterPoints = new HashSet<>();
                newClusterPoints.add(i);
                pointsInClusterList.add(newClusterPoints);
                pointToLeaderMap[i] = newLeaderId;

                Set<Integer> pointsToReEvaluate = new HashSet<>(pointsInClusterList.get(bestLeaderIndex));
                for (int pointIndexToReEvaluate : pointsToReEvaluate) {
                    long[] p = pageData[pointIndexToReEvaluate];
                    long costToOldLeader = calculateTotalLogCost(p, medoidsList.get(bestLeaderIndex), dim);
                    long costToNewLeader = calculateTotalLogCost(p, currentPoint, dim);

                    if (costToNewLeader < costToOldLeader) {
                        pointsInClusterList.get(bestLeaderIndex).remove(pointIndexToReEvaluate);
                        pointsInClusterList.get(newLeaderId).add(pointIndexToReEvaluate);
                        pointToLeaderMap[pointIndexToReEvaluate] = newLeaderId;
                    }
                }
            } else {
                pointsInClusterList.get(bestLeaderIndex).add(i);
                pointToLeaderMap[i] = bestLeaderIndex;
            }
        }

        int k = medoidsList.size();
        long[][] discoveredMedoids = medoidsList.toArray(new long[0][]);
        long[] rawClusterSizes = new long[k];
        for (int i = 0; i < k; i++) {
            rawClusterSizes[i] = pointsInClusterList.get(i).size();
        }



        return sortResults(discoveredMedoids, pointToLeaderMap, rawClusterSizes);
    }

    public static Object[] adaptiveGreedyBasePointSelection_simd(long[][] pageData, int dim) {



        int n = pageData.length;
        if (n == 0) {
            return new Object[]{new long[0][dim], new int[0], new long[0]};
        }


        List<long[]> medoidsList = new ArrayList<>();
        Set<ArrayWrapper> existingMedoidsSet = new HashSet<>();
        List<Set<Integer>> pointsInClusterList = new ArrayList<>();
        int[] pointToLeaderMap = new int[n];

        long[] firstPoint = Arrays.copyOf(pageData[0], dim);
        medoidsList.add(firstPoint);
        existingMedoidsSet.add(new ArrayWrapper(firstPoint));

        Set<Integer> firstClusterPoints = new HashSet<>();
        firstClusterPoints.add(0);
        pointsInClusterList.add(firstClusterPoints);
        pointToLeaderMap[0] = 0;

        for (int i = 1; i < n; i++) {
            long[] currentPoint = pageData[i];

            if (existingMedoidsSet.contains(new ArrayWrapper(currentPoint))) {
                int existingLeaderIndex = -1;
                for (int j = 0; j < medoidsList.size(); j++) {
                    if (Arrays.equals(medoidsList.get(j), currentPoint)) {
                        existingLeaderIndex = j;
                        break;
                    }
                }
                if (existingLeaderIndex != -1) {
                    pointsInClusterList.get(existingLeaderIndex).add(i);
                    pointToLeaderMap[i] = existingLeaderIndex;
                }
                continue;
            }

            int bestLeaderIndex = -1;
            long minCostToExistingLeader = Long.MAX_VALUE;
            for (int j = 0; j < medoidsList.size(); j++) {
                long cost = calculateTotalLogCost_Ultimate(currentPoint, medoidsList.get(j), dim);
                if (cost < minCostToExistingLeader) {
                    minCostToExistingLeader = cost;
                    bestLeaderIndex = j;
                }
            }

            long savingsFromCurrentPoint = minCostToExistingLeader;
            long savingsFromReassignment = 0;

            Set<Integer> pointsInBestLeaderCluster = pointsInClusterList.get(bestLeaderIndex);
            for (int pointIndexInCluster : pointsInBestLeaderCluster) {
                long[] p = pageData[pointIndexInCluster];
                long costToOldLeader = calculateTotalLogCost_Ultimate(p, medoidsList.get(bestLeaderIndex), dim);
                long costToNewPotentialLeader = calculateTotalLogCost_Ultimate(p, currentPoint, dim);
                if (costToNewPotentialLeader < costToOldLeader) {
                    savingsFromReassignment += (costToOldLeader - costToNewPotentialLeader);
                }
            }
            long totalSavings = savingsFromCurrentPoint + savingsFromReassignment;
            long storageCostForNewPoint = calculateBasePointStorageCost_simd(currentPoint, dim);
            if (totalSavings > storageCostForNewPoint) {
                int newLeaderId = medoidsList.size();
                medoidsList.add(Arrays.copyOf(currentPoint, dim));
                existingMedoidsSet.add(new ArrayWrapper(currentPoint));

                Set<Integer> newClusterPoints = new HashSet<>();
                newClusterPoints.add(i);
                pointsInClusterList.add(newClusterPoints);
                pointToLeaderMap[i] = newLeaderId;

                Set<Integer> pointsToReEvaluate = new HashSet<>(pointsInClusterList.get(bestLeaderIndex));
                for (int pointIndexToReEvaluate : pointsToReEvaluate) {
                    long[] p = pageData[pointIndexToReEvaluate];
                    long costToOldLeader = calculateTotalLogCost_Ultimate(p, medoidsList.get(bestLeaderIndex), dim);
                    long costToNewLeader = calculateTotalLogCost_Ultimate(p, currentPoint, dim);

                    if (costToNewLeader < costToOldLeader) {
                        pointsInClusterList.get(bestLeaderIndex).remove(pointIndexToReEvaluate);
                        pointsInClusterList.get(newLeaderId).add(pointIndexToReEvaluate);
                        pointToLeaderMap[pointIndexToReEvaluate] = newLeaderId;
                    }
                }
            } else {
                pointsInClusterList.get(bestLeaderIndex).add(i);
                pointToLeaderMap[i] = bestLeaderIndex;
            }
        }

        int k = medoidsList.size();
        long[][] discoveredMedoids = medoidsList.toArray(new long[0][]);
        long[] rawClusterSizes = new long[k];
        for (int i = 0; i < k; i++) {
            rawClusterSizes[i] = pointsInClusterList.get(i).size();
        }



        return sortResults(discoveredMedoids, pointToLeaderMap, rawClusterSizes);
    }



    record ReassignmentCandidate(int pointIndex, long costToOldLeader, long costToNewPotentialLeader) {}


    public static Object[] adaptiveGreedyBasePointSelection_Ultimate(long[][] pageData, int dim) {

        int n = pageData.length;
        if (n == 0) {
            return new Object[]{new long[0][dim], new int[0], new long[0]};
        }

        List<long[]> medoidsList = new ArrayList<>();
        Set<ArrayWrapper> existingMedoidsSet = new HashSet<>();
        List<Set<Integer>> pointsInClusterList = new ArrayList<>();
        int[] pointToLeaderMap = new int[n];

        
        long[] firstPoint = Arrays.copyOf(pageData[0], dim);
        medoidsList.add(firstPoint);
        existingMedoidsSet.add(new ArrayWrapper(firstPoint));
        Set<Integer> firstClusterPoints = new HashSet<>();
        firstClusterPoints.add(0);
        pointsInClusterList.add(firstClusterPoints);
        pointToLeaderMap[0] = 0;

        
        for (int i = 1; i < n; i++) {
            long[] currentPoint = pageData[i];

            if (existingMedoidsSet.contains(new ArrayWrapper(currentPoint))) {
                int existingLeaderIndex = -1;
                for (int j = 0; j < medoidsList.size(); j++) {
                    if (Arrays.equals(medoidsList.get(j), currentPoint)) {
                        existingLeaderIndex = j;
                        break;
                    }
                }
                if (existingLeaderIndex != -1) {
                    pointsInClusterList.get(existingLeaderIndex).add(i);
                    pointToLeaderMap[i] = existingLeaderIndex;
                }
                continue;
            }

            
            int bestLeaderIndex = -1;
            long minCostToExistingLeader = Long.MAX_VALUE;
            for (int j = 0; j < medoidsList.size(); j++) {
                long cost = calculateTotalLogCost_Ultimate(currentPoint, medoidsList.get(j), dim);
                if (cost < minCostToExistingLeader) {
                    minCostToExistingLeader = cost;
                    bestLeaderIndex = j;
                }
            }

            
            long savingsFromCurrentPoint = minCostToExistingLeader;
            long savingsFromReassignment = 0;

            
            List<ReassignmentCandidate> candidates = new ArrayList<>();
            Set<Integer> pointsInBestLeaderCluster = pointsInClusterList.get(bestLeaderIndex);

            for (int pointIndexInCluster : pointsInBestLeaderCluster) {
                long[] p = pageData[pointIndexInCluster];
                long costToOldLeader = calculateTotalLogCost_Ultimate(p, medoidsList.get(bestLeaderIndex), dim);
                long costToNewPotentialLeader = calculateTotalLogCost_Ultimate(p, currentPoint, dim);

                
                candidates.add(new ReassignmentCandidate(pointIndexInCluster, costToOldLeader, costToNewPotentialLeader));

                if (costToNewPotentialLeader < costToOldLeader) {
                    savingsFromReassignment += (costToOldLeader - costToNewPotentialLeader);
                }
            }
            long totalSavings = savingsFromCurrentPoint + savingsFromReassignment;

            
            long storageCostForNewPoint = calculateBasePointStorageCost_simd(currentPoint, dim);
            if (totalSavings > storageCostForNewPoint) {
                
                int newLeaderId = medoidsList.size();
                medoidsList.add(Arrays.copyOf(currentPoint, dim));
                existingMedoidsSet.add(new ArrayWrapper(currentPoint));

                Set<Integer> newClusterPoints = new HashSet<>();
                newClusterPoints.add(i);
                pointsInClusterList.add(newClusterPoints);
                pointToLeaderMap[i] = newLeaderId;

                
                Set<Integer> oldCluster = pointsInClusterList.get(bestLeaderIndex);
                for (ReassignmentCandidate candidate : candidates) {
                    
                    if (candidate.costToNewPotentialLeader() < candidate.costToOldLeader()) {
                        oldCluster.remove(candidate.pointIndex());               
                        newClusterPoints.add(candidate.pointIndex());           
                        pointToLeaderMap[candidate.pointIndex()] = newLeaderId; 
                    }
                }
            } else {
                
                pointsInClusterList.get(bestLeaderIndex).add(i);
                pointToLeaderMap[i] = bestLeaderIndex;
            }
        }

        
        int k = medoidsList.size();
        long[][] discoveredMedoids = medoidsList.toArray(new long[0][]);
        long[] rawClusterSizes = new long[k];
        for (int i = 0; i < k; i++) {
            rawClusterSizes[i] = pointsInClusterList.get(i).size();
        }

        return sortResults(discoveredMedoids, pointToLeaderMap, rawClusterSizes);
    }

    public static Object[] adaptiveGreedyBasePointSelection_FinalCut(long[][] pageData, int dim) {
        int n = pageData.length;
        if (n == 0) {
            return new Object[]{new long[0][dim], new int[0], new long[0]};
        }

        List<long[]> medoidsList = new ArrayList<>();
        Set<ArrayWrapper> existingMedoidsSet = new HashSet<>();
        List<Set<Integer>> pointsInClusterList = new ArrayList<>();
        int[] pointToLeaderMap = new int[n];

        
        long[] firstPoint = Arrays.copyOf(pageData[0], dim);
        medoidsList.add(firstPoint);
        existingMedoidsSet.add(new ArrayWrapper(firstPoint));
        Set<Integer> firstClusterPoints = new HashSet<>();
        firstClusterPoints.add(0);
        pointsInClusterList.add(firstClusterPoints);
        pointToLeaderMap[0] = 0;

        
        for (int i = 1; i < n; i++) {
            long[] currentPoint = pageData[i];

            if (existingMedoidsSet.contains(new ArrayWrapper(currentPoint))) {
                int existingLeaderIndex = -1;
                for (int j = 0; j < medoidsList.size(); j++) {
                    if (Arrays.equals(medoidsList.get(j), currentPoint)) {
                        existingLeaderIndex = j;
                        break;
                    }
                }
                if (existingLeaderIndex != -1) {
                    pointsInClusterList.get(existingLeaderIndex).add(i);
                    pointToLeaderMap[i] = existingLeaderIndex;
                }
                continue;
            }

            
            int bestLeaderIndex = -1;
            long minCostToExistingLeader = Long.MAX_VALUE;
            for (int j = 0; j < medoidsList.size(); j++) {
                long cost = calculateTotalLogCost_Ultimate(currentPoint, medoidsList.get(j), dim);
                if (cost < minCostToExistingLeader) {
                    minCostToExistingLeader = cost;
                    bestLeaderIndex = j;
                }
            }

            
            long storageCostForNewPoint = calculateBasePointStorageCost_simd(currentPoint, dim);
            long totalSavings = minCostToExistingLeader; 
            boolean worthIt = false;

            
            
            if (totalSavings > storageCostForNewPoint) {
                worthIt = true;
            } else {
                
                Set<Integer> pointsInBestLeaderCluster = pointsInClusterList.get(bestLeaderIndex);
                for (int pointIndexInCluster : pointsInBestLeaderCluster) {
                    long[] p = pageData[pointIndexInCluster];
                    long costToOldLeader = calculateTotalLogCost_Ultimate(p, medoidsList.get(bestLeaderIndex), dim);
                    long costToNewPotentialLeader = calculateTotalLogCost_Ultimate(p, currentPoint, dim);

                    if (costToNewPotentialLeader < costToOldLeader) {
                        totalSavings += (costToOldLeader - costToNewPotentialLeader);
                        
                        if (totalSavings > storageCostForNewPoint) {
                            worthIt = true;
                            break; 
                        }
                    }
                }
            }

            
            if (worthIt) {
                
                int newLeaderId = medoidsList.size();
                medoidsList.add(Arrays.copyOf(currentPoint, dim));
                existingMedoidsSet.add(new ArrayWrapper(currentPoint));

                Set<Integer> newClusterPoints = new HashSet<>();
                newClusterPoints.add(i);
                pointsInClusterList.add(newClusterPoints);
                pointToLeaderMap[i] = newLeaderId;

                
                
                
                Set<Integer> pointsToReEvaluate = new HashSet<>(pointsInClusterList.get(bestLeaderIndex));
                for (int pointIndexToReEvaluate : pointsToReEvaluate) {
                    long[] p = pageData[pointIndexToReEvaluate];
                    long costToOldLeader = calculateTotalLogCost_Ultimate(p, medoidsList.get(bestLeaderIndex), dim);
                    long costToNewLeader = calculateTotalLogCost_Ultimate(p, currentPoint, dim);

                    if (costToNewLeader < costToOldLeader) {
                        pointsInClusterList.get(bestLeaderIndex).remove(pointIndexToReEvaluate);
                        pointsInClusterList.get(newLeaderId).add(pointIndexToReEvaluate);

                        pointToLeaderMap[pointIndexToReEvaluate] = newLeaderId;
                    }
                }
            } else {
                
                pointsInClusterList.get(bestLeaderIndex).add(i);
                pointToLeaderMap[i] = bestLeaderIndex;
            }
        }

        
        int k = medoidsList.size();
        long[][] discoveredMedoids = medoidsList.toArray(new long[0][]);
        long[] rawClusterSizes = new long[k];
        for (int i = 0; i < k; i++) {
            rawClusterSizes[i] = pointsInClusterList.get(i).size();
        }

        
        return sortResults(discoveredMedoids, pointToLeaderMap, rawClusterSizes);
    }




    public static Object[] adaptiveGreedyBasePointSelection_Optimized(long[][] pageData, int dim) {
        int n = pageData.length;
        if (n == 0) {
            return new Object[]{new long[0][dim], new int[0], new long[0]};
        }

        List<long[]> medoidsList = new ArrayList<>();
        Set<ArrayWrapper> existingMedoidsSet = new HashSet<>();
        List<Set<Integer>> pointsInClusterList = new ArrayList<>();
        int[] pointToLeaderMap = new int[n];

        
        medoidsList.add(Arrays.copyOf(pageData[0], dim));
        existingMedoidsSet.add(new ArrayWrapper(pageData[0]));
        Set<Integer> firstClusterPoints = new HashSet<>();
        firstClusterPoints.add(0);
        pointsInClusterList.add(firstClusterPoints);
        pointToLeaderMap[0] = 0;

        
        for (int i = 1; i < n; i++) {
            long[] currentPoint = pageData[i];

            if (existingMedoidsSet.contains(new ArrayWrapper(currentPoint))) {
                
                continue;
            }

            
            int bestLeaderIndex = 0;
            long minCostToExistingLeader = Long.MAX_VALUE;
            for (int j = 0; j < medoidsList.size(); j++) {
                
                long cost = calculateTotalLogCost(currentPoint, medoidsList.get(j), dim);
                if (cost < minCostToExistingLeader) {
                    minCostToExistingLeader = cost;
                    bestLeaderIndex = j;
                }
            }

            
            long savingsFromCurrentPoint = minCostToExistingLeader;
            long savingsFromReassignment = 0;
            Set<Integer> pointsInBestLeaderCluster = pointsInClusterList.get(bestLeaderIndex);
            long[] bestLeaderPoint = medoidsList.get(bestLeaderIndex); 

            for (int pointIndexInCluster : pointsInBestLeaderCluster) {
                long[] p = pageData[pointIndexInCluster];
                
                long costToOldLeader = calculateTotalLogCost(p, bestLeaderPoint, dim);
                long costToNewPotentialLeader = calculateTotalLogCost(p, currentPoint, dim);
                if (costToNewPotentialLeader < costToOldLeader) {
                    savingsFromReassignment += (costToOldLeader - costToNewPotentialLeader);
                }
            }
            long totalSavings = savingsFromCurrentPoint + savingsFromReassignment;

            
            long storageCostForNewPoint = calculateBasePointStorageCost(currentPoint, dim); 
            if (totalSavings > storageCostForNewPoint) {
                
                int newLeaderId = medoidsList.size();
                medoidsList.add(Arrays.copyOf(currentPoint, dim));
                existingMedoidsSet.add(new ArrayWrapper(currentPoint));

                Set<Integer> newClusterPoints = new HashSet<>();
                newClusterPoints.add(i);
                pointsInClusterList.add(newClusterPoints);
                pointToLeaderMap[i] = newLeaderId;

                
                
                Set<Integer> pointsToReEvaluate = new HashSet<>(pointsInClusterList.get(bestLeaderIndex));
                for (int pointIndexToReEvaluate : pointsToReEvaluate) {
                    long[] p = pageData[pointIndexToReEvaluate];
                    
                    long costToNewLeader = calculateTotalLogCost(p, currentPoint, dim);
                    long costToOldLeader = calculateTotalLogCost(p, bestLeaderPoint, dim); 

                    if (costToNewLeader < costToOldLeader) {
                        pointsInClusterList.get(bestLeaderIndex).remove(pointIndexToReEvaluate);
                        pointsInClusterList.get(newLeaderId).add(pointIndexToReEvaluate);
                        pointToLeaderMap[pointIndexToReEvaluate] = newLeaderId;
                    }
                }
            } else {
                
                pointsInClusterList.get(bestLeaderIndex).add(i);
                pointToLeaderMap[i] = bestLeaderIndex;
            }
        }

        
        int k = medoidsList.size();
        long[][] discoveredMedoids = medoidsList.toArray(new long[0][]);
        long[] rawClusterSizes = new long[k];
        for (int i = 0; i < k; i++) {
            rawClusterSizes[i] = pointsInClusterList.get(i).size();
        }
        return sortResults(discoveredMedoids, pointToLeaderMap, rawClusterSizes);
    }






    private static class ClusterInfo implements Comparable<ClusterInfo> {
        long[] medoid;
        List<Integer> pointIndices; 

        ClusterInfo(long[] medoid) {
            this.medoid = medoid;
            this.pointIndices = new ArrayList<>();
        }

        public int getSize() {
            return this.pointIndices.size();
        }

        @Override
        public int compareTo(ClusterInfo other) {
            
            return Integer.compare(this.getSize(), other.getSize());
        }
    }


    public static Object[] adaptiveLeaderBasePointSelectionOpIndex(long[][] pageData, int dim) {
        int n = pageData.length;
        if (n == 0) {
            return new Object[]{new long[0][dim], new int[0], new long[0]};
        }

        
        
        
        
        List<ClusterInfo> clusters = new ArrayList<>();

        
        clusters.add(new ClusterInfo(Arrays.copyOf(pageData[0], dim)));

        for (int i = 0; i < n; i++) {
            long[] currentPoint = pageData[i];
            long minCost = Long.MAX_VALUE;
            int bestLeaderIndex = -1;

            
            for (int j = 0; j < clusters.size(); j++) {
                long cost = calculateTotalLogCost(currentPoint, clusters.get(j).medoid, dim);
                if (cost < minCost) {
                    minCost = cost;
                    bestLeaderIndex = j;
                }
            }

            
            long basePointStorageCost = 0;
            for (int d = 0; d < dim; d++) {
                basePointStorageCost += logEncodingCost(currentPoint[d]) + 1;
            }

            if (minCost > basePointStorageCost) {
                
                ClusterInfo newCluster = new ClusterInfo(Arrays.copyOf(currentPoint, dim));
                newCluster.pointIndices.add(i); 
                clusters.add(newCluster);
            } else {
                
                clusters.get(bestLeaderIndex).pointIndices.add(i);
            }
        }

        
        
        
        
        Collections.sort(clusters);

        
        
        
        
        int k = clusters.size();
        long[][] sortedMedoids = new long[k][dim];
        long[] sortedClusterSize = new long[k];
        int[] sortedClusterAssignment = new int[n];

        for (int newIndex = 0; newIndex < k; newIndex++) {
            ClusterInfo currentCluster = clusters.get(newIndex);

            
            sortedMedoids[newIndex] = currentCluster.medoid;
            sortedClusterSize[newIndex] = currentCluster.getSize();

            
            
            for (int pointIndex : currentCluster.pointIndices) {
                sortedClusterAssignment[pointIndex] = newIndex;
            }
        }

        return new Object[]{sortedMedoids, sortedClusterAssignment, sortedClusterSize};
    }


    public static Object[] adaptiveGreedyBasePointSelection1D(long[][] pageData) {
        long start_t = System.nanoTime();
        int n = pageData.length;
        if (n == 0) {
            return new Object[]{new long[0][1], new int[0], new long[0]};
        }
        int dim = 1;

        class DataPoint implements Comparable<DataPoint> {
            long value;
            int originalIndex;
            DataPoint(long value, int originalIndex) { this.value = value; this.originalIndex = originalIndex; }
            @Override
            public int compareTo(DataPoint other) { return Long.compare(this.value, other.value); }
        }

        DataPoint[] sortedData = new DataPoint[n];
        for (int i = 0; i < n; i++) {
            sortedData[i] = new DataPoint(pageData[i][0], i);
        }
        Arrays.sort(sortedData);

        List<Long> medoids = new ArrayList<>();
        List<List<Integer>> clusterMembers = new ArrayList<>();

        medoids.add(sortedData[0].value);
        clusterMembers.add(new ArrayList<>(Collections.singletonList(0)));

        for (int i = 1; i < n; i++) {
            long currentPointValue = sortedData[i].value;

            int insertionPoint = Collections.binarySearch(medoids, currentPointValue);
            int leaderIdx = -1;
            if (insertionPoint >= 0) { 
                leaderIdx = insertionPoint;
            } else {
                insertionPoint = -(insertionPoint + 1);
                long costBefore = (insertionPoint > 0) ? Math.abs(currentPointValue - medoids.get(insertionPoint - 1)) : Long.MAX_VALUE;
                long costAfter = (insertionPoint < medoids.size()) ? Math.abs(currentPointValue - medoids.get(insertionPoint)) : Long.MAX_VALUE;
                leaderIdx = (costBefore <= costAfter) ? insertionPoint - 1 : insertionPoint;
            }
            long leaderValue = medoids.get(leaderIdx);


            long savingsFromCurrentPoint = logEncodingCost(currentPointValue-leaderValue);


            long savingsFromReassignment = 0;
            for (int memberSortedIndex : clusterMembers.get(leaderIdx)) {
                long memberValue = sortedData[memberSortedIndex].value;
                long costWithOldLeader = logEncodingCost(memberValue-leaderValue);
                long costWithNewLeader = logEncodingCost(memberValue-currentPointValue);
                if (costWithNewLeader < costWithOldLeader) {
                    savingsFromReassignment += (costWithOldLeader - costWithNewLeader);
                }
            }

            long totalSavings = savingsFromCurrentPoint + savingsFromReassignment;
            long storageCostOfNewBasePoint = logEncodingCost(currentPointValue);

            
            if (totalSavings > storageCostOfNewBasePoint) {
                int newMedoidInsertionPoint = Collections.binarySearch(medoids, currentPointValue);
                if (newMedoidInsertionPoint < 0) newMedoidInsertionPoint = -(newMedoidInsertionPoint + 1);

                medoids.add(newMedoidInsertionPoint, currentPointValue);
                clusterMembers.add(newMedoidInsertionPoint, new ArrayList<>());


                List<Integer> pointsToMove = new ArrayList<>();
                int oldLeaderAdjustedIdx = (leaderIdx >= newMedoidInsertionPoint) ? leaderIdx + 1 : leaderIdx;

                for (int memberSortedIndex : clusterMembers.get(oldLeaderAdjustedIdx)) {
                    long memberValue = sortedData[memberSortedIndex].value;
                    if (logEncodingCost(memberValue-currentPointValue)<logEncodingCost(memberValue-leaderValue)) {
                        pointsToMove.add(memberSortedIndex);
                    }
                }

                
                clusterMembers.get(oldLeaderAdjustedIdx).removeAll(pointsToMove);
                clusterMembers.get(newMedoidInsertionPoint).addAll(pointsToMove);

                
                clusterMembers.get(newMedoidInsertionPoint).add(i);

            } else {
                
                clusterMembers.get(leaderIdx).add(i);
            }
        }

        
        int k = medoids.size();
        long[][] discoveredMedoids = new long[k][dim];
        long[] rawClusterSizes = new long[k];
        int[] pointToLeaderMap = new int[n]; 

        for (int i = 0; i < k; i++) {
            discoveredMedoids[i][0] = medoids.get(i);
            rawClusterSizes[i] = clusterMembers.get(i).size();
            for (int memberSortedIndex : clusterMembers.get(i)) {
                pointToLeaderMap[memberSortedIndex] = i;
            }
        }

        
        
        class MedoidInfo {
            long[] medoid; long size; int originalSortedIndex;
            MedoidInfo(long[] medoid, long size, int idx) { this.medoid = medoid; this.size = size; this.originalSortedIndex = idx; }
        }
        List<MedoidInfo> medoidInfoList = new ArrayList<>();
        for(int i=0; i<k; i++) medoidInfoList.add(new MedoidInfo(discoveredMedoids[i], rawClusterSizes[i], i));

        medoidInfoList.sort((a,b) -> Long.compare(b.size, a.size));

        int[] oldToNewIndexMap = new int[k];
        for(int i=0; i<k; i++) oldToNewIndexMap[medoidInfoList.get(i).originalSortedIndex] = i;

        long[][] sortedMedoids = new long[k][dim];
        long[] sortedClusterSizes = new long[k];
        int[] finalClusterAssignment = new int[n];

        for(int i=0; i<k; i++) {
            sortedMedoids[i] = medoidInfoList.get(i).medoid;
            sortedClusterSizes[i] = medoidInfoList.get(i).size;
        }

        for(int i=0; i<n; i++) {
            int originalDataIndex = sortedData[i].originalIndex;
            int leaderInSortedMedoids = pointToLeaderMap[i];
            int finalLeaderIndex = oldToNewIndexMap[leaderInSortedMedoids];
            finalClusterAssignment[originalDataIndex] = finalLeaderIndex;
        }
        System.out.println("time cost: " + (System.nanoTime()-start_t)/8926);
        return new Object[]{sortedMedoids, finalClusterAssignment, sortedClusterSizes};
    }



    private static long[] updateMedoid1D(List<long[]> clusterPoints) {
        
        if (clusterPoints == null || clusterPoints.isEmpty()) {
            return null;
        }
        if (clusterPoints.size() == 1) {
            return clusterPoints.get(0);
        }

        
        clusterPoints.sort(Comparator.comparingLong(p -> p[0]));

        
        int medianIndex = clusterPoints.size() / 2;

        return clusterPoints.get(medianIndex);
    }

    private static long[][] acceleratedInitialization1D(long[][] data, int k) {
        int n = data.length;
        long[][] medoids = new long[k][1];
        Random random = new Random();

        medoids[0] = data[random.nextInt(n)].clone();
        double[] minCosts = new double[n];
        Arrays.fill(minCosts, Double.POSITIVE_INFINITY);

        for (int i = 1; i < k; i++) {
            long newMedoidVal = medoids[i - 1][0];
            double totalCost = 0;

            for (int j = 0; j < n; j++) {
                
                double cost = fastLogEncodingCost(Math.abs(data[j][0] - newMedoidVal));
                minCosts[j] = Math.min(minCosts[j], cost);
                totalCost += minCosts[j];
            }

            double randomVal = random.nextDouble() * totalCost;
            double cumulativeCost = 0;
            for (int j = 0; j < n; j++) {
                cumulativeCost += minCosts[j];
                if (cumulativeCost >= randomVal) {
                    medoids[i] = data[j].clone();
                    break;
                }
            }
        }
        return medoids;
    }

    public static Object[] kMedoidLogCost1D(long[][] data, int k, int maxIter) {
        int n = data.length;
        if (n == 0) return new Object[]{new long[0][1], new int[0], new long[0]};

        
        Set<Long> uniquePoints = new HashSet<>();
        for (long[] point : data) uniquePoints.add(point[0]);
        int distinctCount = uniquePoints.size();
        if (distinctCount < k) {
            k = distinctCount;
        }
        if (k == 0) return new Object[]{new long[0][1], new int[n], new long[0]};

        long[][] medoids = acceleratedInitialization1D(data, k);
        int[] clusterAssignment = new int[n];
        boolean changed = true;

        
        for (int iteration = 0; iteration < maxIter && changed; iteration++) {
            changed = false;

            
            class MedoidWithIndex {
                final long value;
                final int originalIndex;
                MedoidWithIndex(long v, int i) { this.value = v; this.originalIndex = i; }
            }
            List<MedoidWithIndex> sortedMedoidsList = new ArrayList<>(k);
            for (int i = 0; i < k; i++) sortedMedoidsList.add(new MedoidWithIndex(medoids[i][0], i));
            sortedMedoidsList.sort(Comparator.comparingLong(m -> m.value));

            long[] sortedMedoidValues = new long[k];
            int[] sortedOriginalIndices = new int[k];
            for(int i=0; i<k; i++){
                sortedMedoidValues[i] = sortedMedoidsList.get(i).value;
                sortedOriginalIndices[i] = sortedMedoidsList.get(i).originalIndex;
            }

            for (int i = 0; i < n; i++) {
                long point = data[i][0];
                int bestMedoidOriginalIndex;
                int insertionPoint = Arrays.binarySearch(sortedMedoidValues, point);

                if (insertionPoint >= 0) {
                    bestMedoidOriginalIndex = sortedOriginalIndices[insertionPoint];
                } else {
                    insertionPoint = -insertionPoint - 1;
                    if (insertionPoint == 0) {
                        bestMedoidOriginalIndex = sortedOriginalIndices[0];
                    } else if (insertionPoint == k) {
                        bestMedoidOriginalIndex = sortedOriginalIndices[k - 1];
                    } else {
                        
                        long costLeft = fastLogEncodingCost(Math.abs(point - sortedMedoidValues[insertionPoint - 1]));
                        long costRight = fastLogEncodingCost(Math.abs(point - sortedMedoidValues[insertionPoint]));
                        bestMedoidOriginalIndex = (costLeft <= costRight) ? sortedOriginalIndices[insertionPoint - 1] : sortedOriginalIndices[insertionPoint];
                    }
                }
                if (clusterAssignment[i] != bestMedoidOriginalIndex) {
                    changed = true;
                }
                clusterAssignment[i] = bestMedoidOriginalIndex;
            }

            if (!changed) break;

            
            Map<Integer, List<Long>> clusterMap = new HashMap<>();
            for (int i = 0; i < n; i++) {
                clusterMap.computeIfAbsent(clusterAssignment[i], key -> new ArrayList<>()).add(data[i][0]);
            }

            for (int i = 0; i < k; i++) {
                List<Long> pointsInCluster = clusterMap.get(i);
                if (pointsInCluster != null && !pointsInCluster.isEmpty()) {
                    Collections.sort(pointsInCluster);
                    medoids[i][0] = pointsInCluster.get(pointsInCluster.size() / 2);
                }
            }
        }

        
        long[] finalClusterSizes = new long[k];
        for (int assignment : clusterAssignment) if (assignment != -1) finalClusterSizes[assignment]++;

        Integer[] indices = new Integer[k];
        for (int i = 0; i < k; i++) indices[i] = i;
        Arrays.sort(indices, (i1, i2) -> Long.compare(finalClusterSizes[i2], finalClusterSizes[i1]));

        long[][] sortedMedoids = new long[k][1];
        long[] sortedClusterSizes = new long[k];
        int[] oldToNewIndexMap = new int[k];

        for (int i = 0; i < k; i++) {
            int oldIndex = indices[i];
            sortedMedoids[i] = medoids[oldIndex];
            sortedClusterSizes[i] = finalClusterSizes[oldIndex];
            oldToNewIndexMap[oldIndex] = i;
        }

        int[] newClusterAssignment = new int[n];
        for (int i = 0; i < n; i++) {
            newClusterAssignment[i] = oldToNewIndexMap[clusterAssignment[i]];
        }

        return new Object[]{sortedMedoids, newClusterAssignment, sortedClusterSizes};
    }

}
