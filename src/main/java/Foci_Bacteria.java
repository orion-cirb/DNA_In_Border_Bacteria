import FociBacteria_Tools.Tools;
import ij.*;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
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
 * Detect bacteria with OmniPose and foci with DOG
 * @author Orion-CIRB
 */
public class Foci_Bacteria implements PlugIn {
    
    Tools tools = new Tools();
    private String imageDir = "";
    public String outDirResults = "";
    public BufferedWriter distResults;
    public BufferedWriter colocResults;
   
    
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
                     "Nb foci ch1\tNb foci ch2\tFocus channel\tFocus ID\tFocus-bacterium pole distance\t" +
                     "Focus-focus 1 distance\tFocus-focus 2 distance\tFocus-focus 3 distance\tFocus-focus 4 distance\t" +
                     "Focus-focus 5 distance\tFocus-focus 6 distance\n";
            FileWriter fwDistResults = new FileWriter(outDirResults + "distances.xls", false);
            distResults = new BufferedWriter(fwDistResults);
            distResults.write(header);
            distResults.flush();
            header = "Image name\tBacterium ID\tBacterium area (µm2)\tBacterium length (µm)\t" +
                     "Nb foci ch1\tNb foci ch2\tColocalization?\tColocalizing focus ch1 ID\tFocus ch1-bacterium pole distance\t" +
                     "Colocalizing focus ch2 ID\tFocus ch2-bacterium pole distance\n";
            FileWriter fwColocResults = new FileWriter(outDirResults + "colocalization.xls", false);
            colocResults = new BufferedWriter(fwColocResults);
            colocResults.write(header);
            colocResults.flush();
            
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
                ImagePlus bactStack = BF.openImagePlus(options)[indexCh];
                ImagePlus imgBact = tools.doZProjection(bactStack, ZProjector.AVG_METHOD);
                tools.flush_close(bactStack);
                
                // Detect bacteria with Omnipose
                tools.print("- Detecting bacteria -");
                Objects3DIntPopulation bactPop = tools.omniposeDetection(imgBact);
                System.out.println(bactPop.getNbObjects() + " bacteria found");
                
                // Open foci1 channel 1
                indexCh = ArrayUtils.indexOf(channels, chs[1]);
                System.out.println("- Opening foci1 channel " + chs[1] + " -");
                ImagePlus foci1Stack = BF.openImagePlus(options)[indexCh];
                ImagePlus imgFoci1 = tools.doZProjection(foci1Stack, ZProjector.MAX_METHOD);
                tools.flush_close(foci1Stack);
                
                // Detect foci1
                tools.print("- Detecting foci1 -");
                Objects3DIntPopulation foci1Pop = tools.findFoci(imgFoci1);
                System.out.println(foci1Pop.getNbObjects() + " foci1 found");
                tools.fociBactLink(bactPop, foci1Pop);
                System.out.println(foci1Pop.getNbObjects() + " foci1 found in bacteria");
                tools.flush_close(imgFoci1);
                
                // Open foci2 channel
                indexCh = ArrayUtils.indexOf(channels, chs[2]);
                System.out.println("- Opening foci2 channel " + chs[2] + " -");
                ImagePlus foci2Stack = BF.openImagePlus(options)[indexCh];
                ImagePlus imgFoci2 = tools.doZProjection(foci2Stack, ZProjector.MAX_METHOD);
                tools.flush_close(foci2Stack);
                
                // Detect foci
                tools.print("- Detecting foci2 -");
                Objects3DIntPopulation foci2Pop = tools.findFoci(imgFoci2);
                System.out.println(foci2Pop.getNbObjects() + " foci2 found");
                tools.fociBactLink(bactPop, foci2Pop);
                System.out.println(foci2Pop.getNbObjects() + " foci2 found in bacteria");
                tools.flush_close(imgFoci2);
                                
                // Save results
                tools.print("- Saving results -");
                tools.saveResults(bactPop, foci1Pop, foci2Pop, rootName, distResults, colocResults);
                
                // Save images
                tools.drawResults(imgBact, bactPop, foci1Pop, foci2Pop, rootName, outDirResults);
                tools.flush_close(imgBact);
            }
        
            tools.print("--- All done! ---");
            
        }   catch (IOException | FormatException | DependencyException | ServiceException ex) {
            Logger.getLogger(Foci_Bacteria.class.getName()).log(Level.SEVERE, null, ex);
        }  
    }
}    
