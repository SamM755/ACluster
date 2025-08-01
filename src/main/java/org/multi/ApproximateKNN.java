package org.multi;

import java.util.*;

public class ApproximateKNN {
    private double[][] clusterCenters; 
    private long[] clusterSizes;      
    private Random rand = new Random();

    public ApproximateKNN(double[][] clusterCenters, long[] clusterSizes) {
        if (clusterCenters.length != clusterSizes.length) {
            throw new IllegalArgumentException("Cluster centers and sizes must have the same length!");
        }
        this.clusterCenters = clusterCenters;
        this.clusterSizes = clusterSizes;
    }


    public double[][] search_d(double[] query, int k) {
        int dim = query.length; 

        
        PriorityQueue<double[]> minHeap = new PriorityQueue<>(Comparator.comparingDouble(a -> a[dim])); 

        for (int i = 0; i < clusterCenters.length; i++) {
            double distance = euclideanDistance(query, clusterCenters[i]); 
            double[] entry = new double[dim + 1]; 
            System.arraycopy(clusterCenters[i], 0, entry, 0, dim);
            entry[dim] = distance;
            minHeap.offer(entry);
        }

        
        double[][] nearestPoints = new double[k][dim];
        int count = 0;

        while (count < k && !minHeap.isEmpty()) {
            double[] closest = minHeap.poll(); 
            int clusterIndex = findClusterIndex(closest, clusterCenters); 
            long numPointsInCluster = clusterSizes[clusterIndex]; 
            System.out.println(Arrays.toString(clusterCenters[clusterIndex]));
            
            for (int i = 0; i < numPointsInCluster && count < k; i++) {
                System.arraycopy(closest, 0, nearestPoints[count], 0, dim);
                count++;
            }
        }

        return nearestPoints;
    }

    public int[] search(double[] query, int k) {
        int dim = query.length;

        
        PriorityQueue<double[]> minHeap = new PriorityQueue<>(Comparator.comparingDouble(a -> a[dim])); 

        for (int i = 0; i < clusterCenters.length; i++) {
            double distance = euclideanDistance(query, clusterCenters[i]); 
            double[] entry = new double[dim + 1]; 
            System.arraycopy(clusterCenters[i], 0, entry, 0, dim);
            entry[dim] = distance;
            minHeap.offer(entry);
        }

        
        List<Integer> selectedIndices = new ArrayList<>();
        long totalSize = 0;

        while (!minHeap.isEmpty() && totalSize <= k) {
            double[] closest = minHeap.poll(); 
            int clusterIndex = findClusterIndex(closest, clusterCenters); 

            System.out.println(Arrays.toString(clusterCenters[clusterIndex]));
            selectedIndices.add(clusterIndex);
            totalSize += clusterSizes[clusterIndex];
        }

        return selectedIndices.stream().mapToInt(i -> i).toArray();
    }


    private static double euclideanDistance(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }


    private int findClusterIndex(double[] closest, double[][] clusterCenters) {
        int dim = clusterCenters[0].length; 
        for (int i = 0; i < clusterCenters.length; i++) {
            if (isSamePoint(closest, clusterCenters[i], dim)) {
                return i;
            }
        }
        return -1; 
    }

    private boolean isSamePoint(double[] a, double[] b, int dim) {
        for (int i = 0; i < dim; i++) {
            if (Double.compare(a[i], b[i]) != 0) { 
                return false;
            }
        }
        return true;
    }

    public static void main(String[] args) {
        
        double[][] clusterCenters = {
                {20.5, 30.2}, {50.3, 70.8}, {100.1, 120.6}, {200.7, 220.4}, {50.3, 70.8}
        };
        long[] clusterSizes = {100, 200, 50, 10, 10}; 

        ApproximateKNN knn = new ApproximateKNN(clusterCenters, clusterSizes);
        double[] queryPoint = {50.3, 70.8}; 
        int k = 30;

        int[] nearestNeighbors = knn.search(queryPoint, k);
        System.out.println("Approximate KNN Neighbors: " + Arrays.toString(nearestNeighbors));


    }
}
