package org.multi;

import java.util.*;

public class Huff {

    
    public static class HuffmanNode implements Comparable<HuffmanNode> {
        long[] symbol; 
        long freq;
        HuffmanNode left, right;

        public HuffmanNode(long[] symbol, long freq) {
            this.symbol = symbol;
            this.freq = freq;
            this.left = null;
            this.right = null;
        }

        @Override
        public int compareTo(HuffmanNode other) {
            return Long.compare(this.freq, other.freq);  
        }
    }

    
    public static HuffmanNode buildHuffmanTree(long[][] symbols, long[] frequencies) {
        PriorityQueue<HuffmanNode> heap = new PriorityQueue<>();

        
        for (int i = 0; i < symbols.length; i++) {
            heap.add(new HuffmanNode(symbols[i], frequencies[i]));
        }

        
        while (heap.size() > 1) {
            HuffmanNode left = heap.poll();
            HuffmanNode right = heap.poll();
            HuffmanNode merged = new HuffmanNode(null, left.freq + right.freq);
            merged.left = left;
            merged.right = right;
            heap.add(merged);
        }

        return heap.poll();
    }

    
    public static String[] generateHuffmanCodes(long[][] symbols, long[] frequencies) {
        
        HuffmanNode root = buildHuffmanTree(symbols, frequencies);
        Map<String, String> huffmanCodes = new HashMap<>();
        generateCodes(root, "", huffmanCodes);

        
        String[] codes = new String[symbols.length];
        for (int i = 0; i < symbols.length; i++) {
            String key = Arrays.toString(symbols[i]);
            codes[i] = huffmanCodes.get(key);
        }

        return codes;
    }

    
    private static void generateCodes(HuffmanNode node, String code, Map<String, String> huffmanCodes) {
        if (node.symbol != null) {
            huffmanCodes.put(Arrays.toString(node.symbol), code);  
            return;
        }
        if (node.left != null) {
            generateCodes(node.left, code + "0", huffmanCodes);  
        }
        if (node.right != null) {
            generateCodes(node.right, code + "1", huffmanCodes);  
        }
    }

    
    public static Object[] encodeLengthsAndSymbols(String[] huffmanCodes, long[][] symbols) {
        int[] lengths = new int[huffmanCodes.length];
        long[][] sortedSymbols = new long[symbols.length][2];

        
        for (int i = 0; i < huffmanCodes.length; i++) {
            lengths[i] = huffmanCodes[i].length();
            sortedSymbols[i] = symbols[i];
        }

        return new Object[]{lengths, sortedSymbols};
    }

    
    public static void main(String[] args) {
        
        long[][] symbols = {
                {1, 2},
                {3, 4},
                {5, 6},
                {7, 8},
                {9, 10},
                {11, 12}
        };
        long[] frequencies = {5, 9, 12, 13, 16, 45};

        
        String[] huffmanCodes = generateHuffmanCodes(symbols, frequencies);
        System.out.println("Huffman Codes: " + Arrays.toString(huffmanCodes));

        
        Object[] lengthSymbolResult = encodeLengthsAndSymbols(huffmanCodes, symbols);
        int[] lengths = (int[]) lengthSymbolResult[0];
        long[][] sortedSymbols = (long[][]) lengthSymbolResult[1];

        System.out.println("Lengths: " + Arrays.toString(lengths));
        System.out.println("Symbols: " + Arrays.deepToString(sortedSymbols));
    }
}
