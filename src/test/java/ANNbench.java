import org.multi.ClusterDecoder;

import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class ANNbench {
    
    
    

    private static double euclideanDistance(double[] p1, double[] p2) {
        double sum = 0.0;
        for (int i = 0; i < p1.length; i++) {
            double diff = p1[i] - p2[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }

    
    private static class PointWithDistance implements Comparable<PointWithDistance> {
        double[] point;
        double distance;
        PointWithDistance(double[] point, double distance) { this.point = point; this.distance = distance; }
        @Override
        public int compareTo(PointWithDistance other) { return Double.compare(this.distance, other.distance); }
    }

    private static class FlatMedoidInfo implements Comparable<FlatMedoidInfo> {
        final int globalIndex;
        final double distance;
        final long frequency;

        public FlatMedoidInfo(int globalIndex, double distance, long frequency) {
            this.globalIndex = globalIndex;
            this.distance = distance;
            this.frequency = frequency;
        }

        @Override
        public int compareTo(FlatMedoidInfo other) {
            return Double.compare(this.distance, other.distance);
        }
    }



    
    
    
    public static List<double[]> exactKnn(double[][] data, double[] query, int k) {
        if (data == null || data.length == 0 || k <= 0) return new ArrayList<>();
        if (k > data.length) k = data.length;

        List<PointWithDistance> distances = new ArrayList<>(data.length);
        for (double[] point : data) {
            distances.add(new PointWithDistance(point, euclideanDistance(query, point)));
        }
        Collections.sort(distances);
        return distances.stream().limit(k).map(pwd -> pwd.point).collect(Collectors.toList());
    }

    private static class MedoidInfo implements Comparable<MedoidInfo> {
        final int globalIndex;
        final double distance;
        
        final int globalDataStartIndex;
        final int pointCount;

        MedoidInfo(int index, double distance, int globalDataStartIndex, int pointCount) {
            this.globalIndex = index;
            this.distance = distance;
            this.globalDataStartIndex = globalDataStartIndex;
            this.pointCount = pointCount;
        }

        @Override
        public int compareTo(MedoidInfo other) {
            return Double.compare(this.distance, other.distance);
        }
    }

    public static List<double[]> pageAnnSearch(
            double[][] finalDecodedData,
            List<double[][]> medoidsPerPage,
            List<long[]> frequenciesPerPage, 
            double[] query, int k, int m) {

        int pageSize = 10000;
        
        List<double[]> aggregatedPageLevelWinners = new ArrayList<>();
        int numPages = medoidsPerPage.size();

        
        for (int pageIdx = 0; pageIdx < numPages; pageIdx++) {

            
            double[][] pageLocalMedoids = medoidsPerPage.get(pageIdx);
            long[] pageLocalFreqs = frequenciesPerPage.get(pageIdx); 
            if (pageLocalMedoids.length == 0) {
                continue; 
            }

            
            int pageGlobalStartIndex = pageIdx * pageSize;

            
            List<MedoidInfo> pageMedoidInfoList = new ArrayList<>();
            int pointsOffsetWithinPage = 0; 

            for (int medoidLocalIdx = 0; medoidLocalIdx < pageLocalMedoids.length; medoidLocalIdx++) {
                double dist = euclideanDistance(query, pageLocalMedoids[medoidLocalIdx]);
                int pointCount = (int) pageLocalFreqs[medoidLocalIdx];
                
                int globalDataStartIndex = pageGlobalStartIndex + pointsOffsetWithinPage;

                pageMedoidInfoList.add(new MedoidInfo(medoidLocalIdx, dist, globalDataStartIndex, pointCount));

                
                pointsOffsetWithinPage += pointCount;
            }
            Collections.sort(pageMedoidInfoList); 

            
            List<MedoidInfo> promisingMedoidsInPage = new ArrayList<>();
            long pointsCovered = 0;
            for (MedoidInfo medoidInfo : pageMedoidInfoList) {
                promisingMedoidsInPage.add(medoidInfo);
                pointsCovered += medoidInfo.pointCount;
                if (pointsCovered > m) { 
                    break;
                }
            }

            
            List<double[]> pageCandidatePoints = new ArrayList<>();
            for (MedoidInfo promisingMedoid : promisingMedoidsInPage) {
                int startIdx = promisingMedoid.globalDataStartIndex;
                int endIdx = startIdx + promisingMedoid.pointCount;
                
                if (endIdx > finalDecodedData.length) {
                    endIdx = finalDecodedData.length;
                }
                for (int i = startIdx; i < endIdx; i++) {
                    pageCandidatePoints.add(finalDecodedData[i]);
                }
            }

            
            if (!pageCandidatePoints.isEmpty()) {
                double[][] pageCandidateArray = pageCandidatePoints.toArray(new double[0][]);
                
                List<double[]> pageTopK = exactKnn(pageCandidateArray, query, k);
                aggregatedPageLevelWinners.addAll(pageTopK);
            }
        }

        
        if (aggregatedPageLevelWinners.isEmpty()) {
            return new ArrayList<>(); 
        }

        double[][] finalCandidateArray = aggregatedPageLevelWinners.toArray(new double[0][]);
        return exactKnn(finalCandidateArray, query, k);
    }
    
    private static class ArrayWrapper {
        private final double[] data;
        private final int hashCode;

        public ArrayWrapper(double[] data) {
            this.data = data;
            this.hashCode = Arrays.hashCode(data); 
        }

        @Override
        public int hashCode() {
            return hashCode;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            ArrayWrapper that = (ArrayWrapper) o;
            return Arrays.equals(this.data, that.data); 
        }
    }

    public static double calculateRecall(List<double[]> exactKnnResult, List<double[]> annResult) {
        
        if (exactKnnResult == null || exactKnnResult.isEmpty()) {
            return (annResult == null || annResult.isEmpty()) ? 1.0 : 0.0;
        }
        if (annResult == null || annResult.isEmpty()) {
            return 0.0;
        }

        
        
        Map<ArrayWrapper, Long> exactFrequencies = exactKnnResult.stream()
                .map(ArrayWrapper::new)
                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));

        Map<ArrayWrapper, Long> annFrequencies = annResult.stream()
                .map(ArrayWrapper::new)
                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));

        
        long intersectionCount = 0;
        
        for (Map.Entry<ArrayWrapper, Long> entry : exactFrequencies.entrySet()) {
            ArrayWrapper pointWrapper = entry.getKey();
            long exactCount = entry.getValue();

            
            long annCount = annFrequencies.getOrDefault(pointWrapper, 0L);

            
            intersectionCount += Math.min(exactCount, annCount);
        }

        
        return (double) intersectionCount / exactKnnResult.size();
    }


    public static double[] generateNoisyQuery(double[][] data, double noiseStdDev) {
        
        Random random = new Random();

        if (data == null || data.length == 0 || data[0].length == 0) {
            throw new IllegalArgumentException("输入数据集不能为空或无效。");
        }

        
        int randomIndex = random.nextInt(data.length);
        double[] originalPoint = data[randomIndex];

        
        double[] noisyPoint = Arrays.copyOf(originalPoint, originalPoint.length);

        
        for (int i = 0; i < noisyPoint.length; i++) {
            
            
            double noise = random.nextGaussian() * noiseStdDev;
            noisyPoint[i] += noise;
        }

        
        return noisyPoint;
    }

    
    
    
    public static void main(String[] args) throws IOException {
        String fileRoot = "src/test/resources/results/";
        String fileEnd = "_compressed.csv";
        int dataset_indx = 0;
        String[] datasets= {"SSD-bench" ,"profile-income" ,"Blockchain-tr_fixed", "Crop", "gas"};

        int k = 100;
        int m = 5*k;


        // loop for all datasets
        for(int j=0; j<1; j++){
            String targetfile = fileRoot + datasets[j] + fileEnd;
            Object[] results = ClusterDecoder.decode_PageByPage(targetfile);
            double[][] finalDecodedData = (double[][]) results[0];
            List<double[][]> finalMedoidsAsDoubleList = (List<double[][]>) results[1];
            List<long[]> decodedFrequencyList  = (List<long[]>) results[2];

            double[][] flatData = finalDecodedData;

            long knnTime = 0;
            long annTime = 0;
            double total_recall = 0;
            for(int i=0; i<1000; i++){
                double[] query = generateNoisyQuery(finalDecodedData, 30);


                long startTime = System.nanoTime();
                List<double[]> exactResult = exactKnn(flatData, query, k);

                long endTime = System.nanoTime();
                long exact_duration = endTime - startTime;
                knnTime += exact_duration;



                startTime = System.nanoTime();
                List<double[]> pageResult = pageAnnSearch(flatData,finalMedoidsAsDoubleList,decodedFrequencyList,query,k,m);

                endTime = System.nanoTime();
                long page_duration = endTime - startTime;
                annTime += page_duration;

                double rcall = calculateRecall(exactResult, pageResult);
                total_recall += rcall;
            }

            System.out.println("Dataset: " + datasets[j]);
            System.out.println(targetfile);
            System.out.println("KNN avg Time: " + knnTime/1000);
            System.out.println("KNN 1000 Time: " + knnTime);
            System.out.println("ANN avg Time: " +  annTime/1000);
            System.out.println("ANN 1000 Time: " + annTime);
            System.out.println("Avg recall: " + total_recall/1000);
        }










































    }
}
