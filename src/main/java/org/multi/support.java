package org.multi;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class support {

    public static double[][] readDataAsDouble(String filePath) throws IOException {
        
        List<double[]> dataList = new ArrayList<>();

        
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            int lineNumber = 1; 

            
            while ((line = br.readLine()) != null) {
                
                if (line.trim().isEmpty()) {
                    continue;
                }

                
                String[] stringValues = line.split(",");

                
                double[] doubleValues = new double[stringValues.length];

                try {
                    
                    for (int i = 0; i < stringValues.length; i++) {
                        
                        doubleValues[i] = Double.parseDouble(stringValues[i].trim());
                    }
                    
                    dataList.add(doubleValues);

                } catch (NumberFormatException e) {
                    
                    throw new IOException(e);
                }
                lineNumber++;
            }
        }

        
        
        return dataList.toArray(new double[0][0]);
    }

    public static byte[] bitarrayToBytes(BitSet bitset) {
        if (bitset == null || bitset.isEmpty()) {
            return new byte[0];
        }

        BitSet newBitset = new BitSet(bitset.length() + 1);
        newBitset.set(0);

        for (int i = 0; i < bitset.length(); i++) {
            newBitset.set(i + 1, bitset.get(i));
        }

        int length = newBitset.length();
        int paddingLength = (8 - (length % 8)) % 8;
        BitSet paddedBitset = new BitSet(length + paddingLength);
        for (int i = 0; i < length; i++) {
            paddedBitset.set(i + paddingLength, newBitset.get(i));
        }

        byte[] bytes = new byte[(paddedBitset.length() + 7) / 8];
        for (int i = 0; i < paddedBitset.length(); i++) {
            if (paddedBitset.get(i)) {
                bytes[i / 8] |= (byte) (1 << (7 - i % 8));
            }
        }

        return bytes;
    }

    public static int bitsRequiredNoSign(long value) {
        if (value == 0) {
            return 1;
        }
        return (int) Math.ceil(Math.log(Math.abs(value) + 1) / Math.log(2));
    }

    public static String[] readDataAsString(String filePath) throws IOException {
        List<String> data = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                
                
                String[] values = line.split(",");
                StringBuilder processedLine = new StringBuilder();

                for (String value : values) {
                    
                    processedLine.append(value.trim()).append(";"); 
                }

                
                if (processedLine.length() > 0) {
                    processedLine.setLength(processedLine.length() - 1);
                }
                data.add(processedLine.toString());
            }
        }

        
        return data.toArray(new String[0]);
    }

    public static double[][] convertToDouble(String[] data, int dim) {
        if (data == null || data.length == 0) {
            throw new IllegalArgumentException("Input data cannot be null or empty.");
        }

        
        double[][] doubleData = new double[data.length][dim];

        for (int i = 0; i < data.length; i++) {
            
            String[] values = data[i].split(";");

            
            if (values.length != dim) {
                throw new IllegalArgumentException("Dimension mismatch at line " + (i + 1) + ". Expected " + dim + ", but got " + values.length);
            }

            
            for (int j = 0; j < dim; j++) {
                doubleData[i][j] = Double.parseDouble(values[j].trim());
            }
        }

        return doubleData;
    }

    public static void appendToBitstream(BitBuffer bitBuffer, long num, int bits) {
        
        for (int i = bits - 1; i >= 0; i--) {
            boolean bit = (num & (1L << i)) != 0;  
            bitBuffer.appendBit(bit);
        }
    }

    public static void appendToBitstreamStr(BitBuffer bitBuffer, String huffCode, int fixedLength) {
        
        if (huffCode.length() > fixedLength) {
            throw new IllegalArgumentException("Huffman code length exceeds fixed length");
        }

        
        int padding = fixedLength - huffCode.length();
        for (int i = 0; i < padding; i++) {
            bitBuffer.appendBit(false);  
        }

        
        for (int i = 0; i < huffCode.length(); i++) {
            boolean bit = huffCode.charAt(i) == '1';
            bitBuffer.appendBit(bit);
        }
    }

    public static void main(String[] args) {
        BitBuffer buffer = new BitBuffer();

        













        System.out.println(bitsRequiredNoSign(150));
    }
}
