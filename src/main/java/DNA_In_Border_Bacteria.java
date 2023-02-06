import DNA_In_Border_Bacteria_Tools.Tools;

import ij.*;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;


/**
 * Detect bacteria with OmniPose
 * Measure DNA intensity inside and on the edges of bacteria
 * 
 * @author ORION-CIRB
 */
public class DNA_In_Border_Bacteria implements PlugIn {
    
    Tools tools = new Tools();
    private String imageDir = "";
    public String outDirResults = "";
    public BufferedWriter results;
   
    
    public void run(String arg) {
        try {
            if (!tools.checkInstalledModules()) {
                return;
            } 
            
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }   
            
            // Find images with extension
            String file_ext = tools.findImageType(new File(imageDir));
            ArrayList<String> imageFiles = tools.findImages(imageDir, file_ext);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("Error", "No images found with " + file_ext + " extension");
                return;
            }
            
            // Create output folder
            outDirResults = imageDir + File.separator + "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            String header = "Image name\tBacterium ID\tBacterium area (µm2)\tBacterium length (µm)\t" +
                     "DNA mean intensity inside bacterium\tDNA mean intensity on bacterium edges\n";
            FileWriter fwDistResults = new FileWriter(outDirResults + "results.xls", false);
            results = new BufferedWriter(fwDistResults);
            results.write(header);
            results.flush();
           
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);
            
            // Find channel names
            String[] channels = tools.findChannels(imageFiles.get(0), meta, reader);

            // Dialog box
            String[] chs = tools.dialog(channels);
            if (chs == null) {
                IJ.showMessage("Error", "Plugin canceled");
                return;
            }
            
            for (String f : imageFiles) {
                reader.setId(f);
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ------");
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setSplitChannels(true);
                
                // Open bacteria channel
                int indexCh = ArrayUtils.indexOf(channels, chs[0]);
                System.out.println("- Opening bacteria channel " + chs[0] + " -");
                ImagePlus stackBact = BF.openImagePlus(options)[indexCh];
                ImagePlus imgBact = (stackBact.getNSlices() == 1) ? new Duplicator().run(stackBact) : new Duplicator().run(stackBact, stackBact.getNSlices()/2, stackBact.getNSlices()/2);
                tools.flush_close(stackBact);
                
                // Detect bacteria with Omnipose
                tools.print("- Detecting bacteria -");
                Objects3DIntPopulation bactPop = tools.omniposeDetection(imgBact);
                System.out.println(bactPop.getNbObjects() + " bacteria found");
                
                // Open DNA channel
                indexCh = ArrayUtils.indexOf(channels, chs[1]);
                System.out.println("- Opening DNA channel " + chs[1] + " -");
                ImagePlus stackDna = BF.openImagePlus(options)[indexCh];
                ImagePlus imgDna = (stackDna.getNSlices() == 1) ? new Duplicator().run(stackDna) : new Duplicator().run(stackDna, stackDna.getNSlices()/2, stackDna.getNSlices()/2);
                tools.flush_close(stackDna);
                
                // Save results
                tools.print("- Saving results -");
                Objects3DIntPopulation bactBorderPop = tools.saveResults(bactPop, imgDna, rootName, outDirResults, results);
                
                // Save images
                tools.drawResults(imgBact, bactPop, "_bacteria.tif", rootName, outDirResults);
                tools.drawResults(imgDna, bactBorderPop, "_edges.tif", rootName, outDirResults);
                
                tools.flush_close(imgBact);
                tools.flush_close(imgDna);
            }
            results.close();
            tools.print("--- All done! ---");
            
        }   catch (IOException | FormatException | DependencyException | ServiceException ex) {
            Logger.getLogger(DNA_In_Border_Bacteria.class.getName()).log(Level.SEVERE, null, ex);
        }  
    }
}    
