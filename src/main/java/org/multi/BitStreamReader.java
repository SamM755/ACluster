package org.multi;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;

public class BitStreamReader {
    private FileInputStream fis;
    private byte[] buffer; 
    private int bufferSize = 0; 
    private int currentBitPos = 0; 
    private int totalBitsInBuffer = 0; 

    
    public BitStreamReader(String filePath) throws IOException {
        this.fis = new FileInputStream(filePath);
        this.buffer = new byte[1024]; 
        refillBuffer();
    }

    
    public byte[] readBytes(int numBytes) throws IOException {
        byte[] result = new byte[numBytes];
        int bytesRead = 0;

        while (bytesRead < numBytes) {
            if (currentBitPos >= totalBitsInBuffer) {
                refillBuffer();
                if (bufferSize == 0) {
                    throw new IOException("Reached end of file while reading bytes.");
                }
            }

            int byteOffset = currentBitPos / 8;
            int bytesToCopy = Math.min(numBytes - bytesRead, bufferSize - byteOffset);
            System.arraycopy(buffer, byteOffset, result, bytesRead, bytesToCopy);

            bytesRead += bytesToCopy;
            currentBitPos += bytesToCopy * 8;
        }

        return result;
    }

    public void skipBits(long numBitsToSkip) throws IOException {
        if (numBitsToSkip < 0) return;

        
        long bitsLeftInBuffer = totalBitsInBuffer - currentBitPos;

        if (numBitsToSkip <= bitsLeftInBuffer) {
            
            currentBitPos += numBitsToSkip;
            return;
        }

        
        
        numBitsToSkip -= bitsLeftInBuffer;

        
        long bytesToSkip = numBitsToSkip / 8;
        if (bytesToSkip > 0) {
            
            
            long actuallySkipped = fis.skip(bytesToSkip);
            if (actuallySkipped != bytesToSkip) {
                
            }
        }

        
        refillBuffer();
        if (bufferSize == 0 && numBitsToSkip > 0) {
            throw new IOException("Reached end of file while trying to skip bits.");
        }

        
        long remainingBitsToSkip = numBitsToSkip % 8;
        if(remainingBitsToSkip > totalBitsInBuffer){
            throw new IOException("Reached end of file while trying to skip remaining bits.");
        }
        currentBitPos = (int) remainingBitsToSkip;
    }

    
    private void refillBuffer() throws IOException {
        bufferSize = fis.read(buffer);
        if (bufferSize > 0) {
            totalBitsInBuffer = bufferSize * 8;
            currentBitPos = 0;
        } else {
            totalBitsInBuffer = 0; 
        }
    }

    
    public void close() throws IOException {
        if (fis != null) {
            fis.close();
        }
    }

    public long readBits(int numBits) throws IOException {
        if (numBits > 64) {
            throw new IllegalArgumentException("Cannot read more than 64 bits at a time.");
        }

        long value = 0;
        for (int i = 0; i < numBits; i++) {
            if (currentBitPos >= totalBitsInBuffer) {
                refillBuffer();
                if (totalBitsInBuffer == 0) {
                    throw new IOException("Reached end of file while reading bits.");
                }
            }

            int byteIndex = currentBitPos / 8;
            int bitIndex = 7 - (currentBitPos % 8);
            boolean bit = (buffer[byteIndex] & (1 << bitIndex)) != 0;

            value = (value << 1) | (bit ? 1L : 0L);
            currentBitPos++;
        }

        return value;
    }

    
    public static long extractBits(byte[] inputBytes, int numBits) {
        if (numBits > 64) {
            throw new IllegalArgumentException("Cannot extract more than 64 bits at a time.");
        }

        long value = 0;
        int totalBits = inputBytes.length * 8;
        int bitPos = 0;

        while (bitPos < numBits && bitPos < totalBits) {
            int byteIndex = bitPos / 8;
            int bitIndex = 7 - (bitPos % 8);
            boolean bit = (inputBytes[byteIndex] & (1 << bitIndex)) != 0;

            value = (value << 1) | (bit ? 1L : 0L);
            bitPos++;
        }

        return value;
    }

    
    public static void main(String[] args) {
        try {
            BitStreamReader reader = new BitStreamReader("C:\\Users\\kenny\\Desktop\\cluster-encode\\cluster-encode\\SSD-bench_compressed_varied.csv");

            
            long first11Bits = reader.readBits(8);
            System.out.println("First 11 Bits: " + Long.toBinaryString(first11Bits));

            
            long next5Bits = reader.readBits(8);
            System.out.println("Next 5 Bits: " + Long.toBinaryString(next5Bits));

            
            byte[] readData = reader.readBytes(2);
            long extractedBits = extractBits(readData, 16);
            System.out.println("Extracted 11 Bits: " + Long.toBinaryString(extractedBits));

            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
