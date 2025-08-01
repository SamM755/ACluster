package org.multi;

import org.apache.poi.ss.usermodel.*;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit; 

public class ExcelWriter {

    public static void writeResultsToExcel(List<Object[]> resultsData, String filePath) {
        
        String[] headers = {"Dataset", "Original Size (Bytes)", "Compressed Size (Bytes)", "Compression Ratio", "Compression Time (ms)"};

        
        try (Workbook workbook = new XSSFWorkbook()) {
            Sheet sheet = workbook.createSheet("Adaptive Compression Results");

            
            CellStyle ratioStyle = workbook.createCellStyle();
            DataFormat dataFormat = workbook.createDataFormat();
            ratioStyle.setDataFormat(dataFormat.getFormat("0.000"));

            
            Row headerRow = sheet.createRow(0);
            for (int i = 0; i < headers.length; i++) {
                Cell cell = headerRow.createCell(i);
                cell.setCellValue(headers[i]);
            }

            
            int rowNum = 1;
            for (Object[] rowData : resultsData) {
                Row row = sheet.createRow(rowNum++);

                row.createCell(0).setCellValue((String) rowData[0]);      
                row.createCell(1).setCellValue((Integer) rowData[1]);     
                row.createCell(2).setCellValue((Long) rowData[2]);        

                Cell ratioCell = row.createCell(3);
                ratioCell.setCellValue((Double) rowData[3]);              
                ratioCell.setCellStyle(ratioStyle); 

                row.createCell(4).setCellValue((Long) rowData[4]);        
            }

            
            for(int i = 0; i < headers.length; i++) {
                sheet.autoSizeColumn(i);
            }

            
            try (FileOutputStream fileOut = new FileOutputStream(filePath)) {
                workbook.write(fileOut);
            }

            System.out.println(filePath);

        } catch (IOException e) {
            System.err.println(e.getMessage());
            e.printStackTrace();
        }
    }
}