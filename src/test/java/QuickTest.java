import org.multi.PageEncoder;
import org.multi.ClusterDecoder;
import java.io.IOException;

public class QuickTest {
    public static void main(String[] args) throws IOException {
        // needs java incubator for java simd
        // run with Java version > 21.0
        // change the run config by adding "--add-modules jdk.incubator.vector" to VM options

        // method_idx
        // for KCluster:0;
        // ACluster:1,
        // below for 1D data only:
        // KOptimal:2
        // AOptimal:3
        int method_idx = 1;

        // change to absolute file path root if needed
        String outputPathRoot = "src/test/resources/results/";
        String fileRoot = "src/test/resources/datasets/";
        int dataset_indx = 4;

        // datafiles needs to be in .csv, no header, columns separated by ","
        String[] datasets= {"SSD-bench" ,"profile-income" ,"Blockchain-tr", "Crop", "gas"};
        String dataset = datasets[dataset_indx];
        int[] dims = {1,1,1,46,64};
        int dim = dims[dataset_indx];
        int[] data_size = {8926,14825,99999, 24000,13910};

        // default bitpacking setups (blocks for reference point, pack for residual)
        int block_size = 10;
        int packSize = 10;

        // default parameters
        int k=1000;
        int pageSize = 10000;

        // for page exp
        int[] pageSizeList= {1000,2000,5000,10000,20000};

        // for K exp
        int[] ks={10,50,100,500,1000};

        // for ann K exp (have to compress accordingly first, then find results using the ANNbench)
        int ann_k[] = {100,500,1000};

        String method_name[] = {"KCluster", "ACluster", "KOptimal", "AOptimal"};
        String filePath = fileRoot + dataset;
        String outputPath = outputPathRoot + dataset;

        int numRows = data_size[dataset_indx];


//        int numByte = data_size[dataset_indx];
        int iterations = 100;
        long final_byte_res = 0;
        long final_time_res = 0;
        long startTime;
        long endTime;
        System.out.println("Encoding...");

        for(int i=0; i<iterations; i++){
            System.out.println("Method name: " + method_name[method_idx]);
            startTime = System.nanoTime();
            Object[] result = PageEncoder.encode(filePath, outputPath, pageSize, k, numRows, dim, packSize, block_size, method_idx);
            endTime = System.nanoTime();
            long duration = endTime - startTime;
            final_byte_res += (long) result[0];
            final_time_res += duration;
//            System.out.println("Time Cost: " + (duration / data_size[dataset_indx]) + " ns per point");
//            System.out.println("Final Byte Cost: " + result[0]);
//            System.out.println("Header bit Cost: " + result[1]);
//            System.out.println("Reference bit Cost: " + result[2]);
//            System.out.println("Reference weights bit Cost: " + result[3]);
//            System.out.println("Residual bit cost: " + result[4]);
        }
        final_byte_res /= iterations;
        final_time_res /= iterations;
        System.out.println("Encode Time: " + (final_time_res / data_size[dataset_indx]) + " ns per point");
        System.out.println("Final Byte: " + final_byte_res);

        long final_decode_time = 0;
        System.out.println("Decoding...");
        for(int i=0; i<iterations; i++){
            startTime = System.nanoTime();
            Object[] results = ClusterDecoder.decode_PageByPage("src/test/resources/" + dataset + "_compressed.csv");
            endTime = System.nanoTime();
            final_decode_time += endTime-startTime;
        }
        final_decode_time /= iterations;
        System.out.println("Decode Time: " + (final_decode_time/data_size[dataset_indx]) + "ns per point");

    }
}
