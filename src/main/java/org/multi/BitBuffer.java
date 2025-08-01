package org.multi;


public class BitBuffer {
    private byte[] buffer;
    private int currentBitPos;

    public BitBuffer() {
        buffer = new byte[1];  
        currentBitPos = 0;     
    }

    public void appendBit(boolean bit) {
        if (currentBitPos >= buffer.length * 8) {
            expandBuffer();
        }

        int byteIndex = currentBitPos / 8;
        int bitIndex = currentBitPos % 8;

        if (bit) {
            buffer[byteIndex] |= (1 << (7 - bitIndex));  
        } else {
            buffer[byteIndex] &= ~(1 << (7 - bitIndex));  
        }

        currentBitPos++;
    }

    private void expandBuffer() {
        byte[] newBuffer = new byte[buffer.length * 2];  
        System.arraycopy(buffer, 0, newBuffer, 0, buffer.length);
        buffer = newBuffer;
    }

    public String toBinaryString() {
        StringBuilder binaryString = new StringBuilder();

        for (int i = 0; i < buffer.length; i++) {
            byte b = buffer[i];
            for (int j = 7; j >= 0; j--) {
                binaryString.append((b >> j) & 1);
            }
        }

        return binaryString.toString().substring(0, currentBitPos);  
    }


public byte[] toByteArray() {
    int bitLength = currentBitPos;

    int paddedBitLength = bitLength;
    while (paddedBitLength % 8 != 0) {
        paddedBitLength++;
    }

    int finalByteLength = paddedBitLength / 8;
    byte[] result = new byte[finalByteLength];

    for (int i = 0; i < bitLength; i++) {
        int originalByteIndex = i / 8;
        int originalBitIndex = i % 8;
        boolean bit = (buffer[originalByteIndex] & (1 << (7 - originalBitIndex))) != 0;

        int newByteIndex = i / 8;
        int newBitIndex = i % 8;

        if (bit) {
            result[newByteIndex] |= (1 << (7 - newBitIndex));
        }
    }

    return result;
}

    public void loadFromByteArray(byte[] byteArray) {
        int totalBits = byteArray.length * 8;

        buffer = new byte[byteArray.length];
        currentBitPos = 0;

        for (int bitPos = 0; bitPos < totalBits; bitPos++) {
            int byteIndex = bitPos / 8;
            int bitIndex = 7 - (bitPos % 8);
            boolean bit = (byteArray[byteIndex] & (1 << bitIndex)) != 0;
            appendBit(bit);
        }
    }



    public void merge(BitBuffer other) {
        for (int i = 0; i < other.currentBitPos; i++) {
            boolean bit = (other.buffer[i / 8] & (1 << (7 - i % 8))) != 0;
            appendBit(bit);
        }
    }

    public int size() {
        return currentBitPos;
    }

    public boolean getBit(int bitIndex) {
        if (bitIndex < 0 || bitIndex >= currentBitPos) {
            throw new IndexOutOfBoundsException("Bit index out of range: " + bitIndex);
        }

        int byteIndex = bitIndex / 8;
        int bitOffset = 7 - (bitIndex % 8);

        return (buffer[byteIndex] & (1 << bitOffset)) != 0;
    }

    public static void main(String[] args) {
        BitBuffer buffer = new BitBuffer();

        String bitString = "11";
        for (char c : bitString.toCharArray()) {
            buffer.appendBit(c == '1');
        }
        System.out.println(buffer.currentBitPos);
        for (byte b : buffer.buffer) {
            System.out.println(String.format("%8s", Integer.toBinaryString(b & 0xFF)).replace(' ', '0'));
        }

        byte[] result = buffer.toByteArray();

        System.out.println("Final Byte Array (Binary):");
        for (byte b : result) {
            System.out.println(String.format("%8s", Integer.toBinaryString(b & 0xFF)).replace(' ', '0'));
        }

        BitBuffer decodedBuffer = new BitBuffer();
        decodedBuffer.loadFromByteArray(result);

        
        System.out.println("Decoded Byte Array (Binary):");
        for (byte b : decodedBuffer.buffer) {
            System.out.println(String.format("%8s", Integer.toBinaryString(b & 0xFF)).replace(' ', '0'));
        }
    }
}

