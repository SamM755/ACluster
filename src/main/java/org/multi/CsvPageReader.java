package org.multi;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class CsvPageReader implements AutoCloseable {

    private final BufferedReader reader;
    private final int pageSize;
    private boolean isFinished = false;


    public CsvPageReader(String filePath, int pageSize) throws IOException {
        if (pageSize <= 0) {
            throw new IllegalArgumentException("Page size must be positive.");
        }
        this.reader = new BufferedReader(new FileReader(filePath));
        this.pageSize = pageSize;
    }


    public double[][] readNextPage() throws IOException {
        if (isFinished) {
            return null; 
        }

        List<double[]> pageDataList = new ArrayList<>(pageSize);
        String line = null;
        int linesRead = 0;

        
        while (linesRead < pageSize && (line = reader.readLine()) != null) {
            if (line.trim().isEmpty()) {
                continue; 
            }

            String[] stringValues = line.split(",");
            double[] doubleValues = new double[stringValues.length];

            try {
                for (int i = 0; i < stringValues.length; i++) {
                    doubleValues[i] = Double.parseDouble(stringValues[i].trim());
                }
                pageDataList.add(doubleValues);
                linesRead++;
            } catch (NumberFormatException e) {
                
                this.close();
                throw new IOException(e);
            }
        }

        
        if (line == null) {
            isFinished = true;
        }

        
        if (pageDataList.isEmpty()) {
            return null;
        }

        return pageDataList.toArray(new double[0][0]);
    }

    @Override
    public void close() throws IOException {
        if (reader != null) {
            reader.close();
        }
    }
}
