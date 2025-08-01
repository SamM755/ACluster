package org.multi;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class PageEncoder {

    public static Object[] encode(String filePath, String outputPathRoot, int pageSize, int k, int data_size, int dim, int packSize, int block_size, int method_idx) throws IOException {
        String data_file = filePath + ".csv";

        String outputFilePath = outputPathRoot + "_compressed.csv";

        System.out.println(outputFilePath);
        int pageCount = 0;
        
        List<int[]> maxDecimalList = new ArrayList<>();
        List<long[]> minValueList = new ArrayList<>();
        List<BitBuffer> frebitList = new ArrayList<>();
        List<BitBuffer> residualbitList = new ArrayList<>();
        List<int[][]> packResMetadataList = new ArrayList<>();
        List<long[]> minBaseList = new ArrayList<>();
        List<int[]> maxBaseBitLenList = new ArrayList<>();
        List<BitBuffer> baseBitstreamList = new ArrayList<>();
        List<Integer> pageDataPointsList = new ArrayList<>();
        List<Integer> pagePackNumList = new ArrayList<>();
        List<Integer> pageKList = new ArrayList<>();


        try (CsvPageReader reader = new CsvPageReader(data_file, pageSize)) {
            double[][] currentPage;
            
            
            while ((currentPage = reader.readNextPage()) != null) {


                
                
                Object[] result = Encoder.clusterEncoderPageNoHuff(currentPage,k,packSize,dim, block_size, filePath, method_idx);
                int[] maxDecimalPlacements = (int[]) result[0];
                long[] minValues = (long[]) result[1];
                BitBuffer frebit = (BitBuffer) result[2];
                BitBuffer resbit = (BitBuffer) result[3];
                int[][] packResMetadata = (int[][]) result[4];
                int pagePackNum = (int) result[5];
                int pageDataPoints = (int) result[6];
                long[] minBase = (long[]) result[7];
                int[] maxBaseBitLen = (int[]) result[8];
                BitBuffer baseBitstream = (BitBuffer) result[9];
                int page_k = (int) result[10];

                maxDecimalList.add(maxDecimalPlacements);
                minValueList.add(minValues);
                pageDataPointsList.add(pageDataPoints);
                pagePackNumList.add(pagePackNum);
                frebitList.add(frebit);
                residualbitList.add(resbit);
                packResMetadataList.add(packResMetadata);
                minBaseList.add(minBase);
                maxBaseBitLenList.add(maxBaseBitLen);
                baseBitstreamList.add(baseBitstream);
                pageKList.add(page_k);

                pageCount++;
            }

        } catch (IOException e) {
            System.err.println(e.getMessage());
            e.printStackTrace();
        }

        return Encoder.saveCompressedDataNoHuff_Page(outputFilePath, dim, packSize, block_size, pageCount,
                maxDecimalList, minValueList, pageDataPointsList, frebitList, residualbitList, packResMetadataList,
                minBaseList, maxBaseBitLenList, baseBitstreamList, pageKList, pageSize);
    }


}
