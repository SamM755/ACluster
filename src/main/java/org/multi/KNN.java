package org.multi;

import java.util.*;

public class KNN {

    public static double[][] knnSearch(double[][] dataset, double[] query, int k) {
        if (dataset == null || dataset.length == 0 || dataset[0].length != query.length || k <= 0) {
            throw new IllegalArgumentException("Dataset cannot be null, empty, or have incorrect dimensions.");
        }

        int dim = query.length; 

        
        PriorityQueue<double[]> maxHeap = new PriorityQueue<>(Comparator.comparingDouble(a -> -a[dim])); 

        for (double[] dataPoint : dataset) {
            double distance = euclideanDistance(query, dataPoint); 

            if (maxHeap.size() < k) {
                double[] entry = new double[dim + 1]; 
                System.arraycopy(dataPoint, 0, entry, 0, dim);
                entry[dim] = distance;
                maxHeap.offer(entry);
            } else if (distance < maxHeap.peek()[dim]) { 
                maxHeap.poll();
                double[] entry = new double[dim + 1];
                System.arraycopy(dataPoint, 0, entry, 0, dim);
                entry[dim] = distance;
                maxHeap.offer(entry);
            }
        }

        
        double[][] nearestNeighbors = new double[k][dim];
        for (int i = k - 1; i >= 0; i--) {
            double[] nearest = Objects.requireNonNull(maxHeap.poll());
            System.arraycopy(nearest, 0, nearestNeighbors[i], 0, dim); 
        }
        


        return nearestNeighbors;
    }

    public static double[][] knnSearch(List<double[]> dataset, double[] query, int k) {
        return knnSearch(dataset.toArray(new double[0][]), query, k);
    }

    private static double euclideanDistance(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }

    public static void main(String[] args) {
        












    }
}
