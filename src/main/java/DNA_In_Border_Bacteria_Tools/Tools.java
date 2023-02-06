package DNA_In_Border_Bacteria_Tools;

import DNA_In_Border_Bacteria.Cellpose.CellposeTaskSettings;
import DNA_In_Border_Bacteria.Cellpose.CellposeSegmentImgPlusAdvanced;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import fiji.util.gui.GenericDialogPlus;
import ij.plugin.RGBStackMerge;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.VoxelInt;
import mcib3d.geom2.measurements.MeasureFeret;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;


/**
 * @author Orion-CIRB
 */
public class Tools {
    private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
      
    public Calibration cal = new Calibration();
    private double pixelSurf = 0;
    String[] channelsName = {"Bacteria: ", "DNA: "};
    
     // Omnipose
    private String omniposeEnvDirPath = IJ.isWindows()? System.getProperty("user.home")+"\\miniconda3\\envs\\omnipose" : "/opt/miniconda3/envs/omnipose";
    private String omniposeModelsPath = IJ.isWindows()? System.getProperty("user.home")+"\\.cellpose\\models\\": System.getProperty("user.home")+"/.cellpose/models/";
    private String omniposeModel = "bact_phase_omnitorch_0";
    private int omniposeDiameter = 0;
    private int omniposeMaskThreshold = 0;
    private double omniposeFlowThreshold = 0;
    private boolean useGpu = true;
    
    // Bacteria
    private double minBactSurface = 1;
    private double maxBactSurface = 20;
    private float bactErosion = 0.4f;
    
    
    
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.showMessage("Error", "3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Flush and close an image
     */
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Find images extension
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
               case "nd" :
                   ext = fileExt;
                   break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "lsm" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
                case "tiff" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }
         
    
    /**
     * Find images in folder
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }

    
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public void findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
    }
    
    
    /**
     * Find channels name
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n);
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;   
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);         
    }
    
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] channels) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String ch : channelsName) {
            gd.addChoice(ch, channels, channels[index]);
            index++;
        }
        
        gd.addMessage("Bacteria detection", Font.getFont("Monospace"), Color.blue);
        gd.addDirectoryField("Omnipose environment directory: ", omniposeEnvDirPath);
        gd.addDirectoryField("Omnipose models path: ", omniposeModelsPath);
        gd.addNumericField("Min bacterium area (µm2): ", minBactSurface);
        gd.addNumericField("Max bacterium area (µm2): ", maxBactSurface);
        gd.addNumericField("Bacterium erosion (µm): ", bactErosion);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm):", cal.pixelWidth);
        gd.showDialog();
        
        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();

        omniposeEnvDirPath = gd.getNextString();
        omniposeModelsPath = gd.getNextString();
        minBactSurface = (float) gd.getNextNumber();
        maxBactSurface = (float) gd.getNextNumber();
        bactErosion = (float) gd.getNextNumber();
        
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = 1;
        pixelSurf = cal.pixelWidth*cal.pixelHeight;
        
        if (gd.wasCanceled())
           ch = null;
                
        return(ch);
    }
    
    
    /**
    * Detect bacteria with Omnipose
    */
    public Objects3DIntPopulation omniposeDetection(ImagePlus imgBact){
        ImagePlus imgIn = new Duplicator().run(imgBact);
        
        // Set Omnipose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(omniposeModelsPath+omniposeModel, 1, omniposeDiameter, omniposeEnvDirPath);
        settings.setVersion("0.7");
        settings.setCluster(true);
        settings.setOmni(true);
        settings.useMxNet(false);
        settings.setCellProbTh(omniposeMaskThreshold);
        settings.setFlowTh(omniposeFlowThreshold);
        settings.useGpu(useGpu);
        
        // Run Omnipose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgIn);
        //PrintStream console = System.out;
        //System.setOut(new NullPrintStream());
        ImagePlus imgOut = cellpose.run();
        //System.setOut(console);
        imgOut.setCalibration(cal);
        
        // Filter bacteria population
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        pop = new Objects3DIntPopulationComputation(pop).getExcludeBorders(ImageHandler.wrap(imgOut), false);
        pop = new Objects3DIntPopulationComputation(pop).getFilterSize(minBactSurface/pixelSurf, maxBactSurface/pixelSurf);
        pop.resetLabels();
        
        // Close images
        flush_close(imgIn);
        flush_close(imgOut);
        
        return(pop);
    }
    

    /**
     * Compute bacteria parameters and save them in file
     * @throws java.io.IOException
     */
    public Objects3DIntPopulation saveResults(Objects3DIntPopulation bactPop, ImagePlus img, String imgName, String outDir, BufferedWriter resFile) throws IOException {
        Objects3DIntPopulation bactBorderPop = new Objects3DIntPopulation();
        for (Object3DInt bact : bactPop.getObjects3DInt()) {
            double bactSurf = new MeasureVolume(bact).getValueMeasurement(MeasureVolume.VOLUME_UNIT);
            VoxelInt feret1Unit = new MeasureFeret(bact).getFeret1Unit();
            VoxelInt feret2Unit = new MeasureFeret(bact).getFeret2Unit();
            double bactLength = feret1Unit.distance(feret2Unit)*cal.pixelWidth;
            
            float erosion = (float)(bactErosion/cal.pixelWidth);
            Object3DInt bactBorder = new Object3DComputation(bact).getObjectEdgeMorpho(erosion, erosion, erosion, false);
            Object3DInt bactInside = new Object3DComputation(bact).getObjectSubtracted(bactBorder);
            
            double volbactInside = new MeasureVolume(bactInside).getValueMeasurement(MeasureVolume.VOLUME_UNIT);
            if (volbactInside != 0) {
                bactBorderPop.addObject(bactBorder);
                double bactInsideInt = new MeasureIntensity(bactInside, ImageHandler.wrap(img)).getValueMeasurement(MeasureIntensity.INTENSITY_AVG);
                double bactBorderInt = new MeasureIntensity(bactBorder, ImageHandler.wrap(img)).getValueMeasurement(MeasureIntensity.INTENSITY_AVG);
                resFile.write(imgName+"\t"+bact.getLabel()+"\t"+bactSurf+"\t"+bactLength+"\t"+bactInsideInt+"\t"+bactBorderInt+"\n");
            } else {
                resFile.write(imgName+"\t"+bact.getLabel()+"\t"+bactSurf+"\t"+bactLength+"\n");
            }
            resFile.flush();
        }
        return bactBorderPop;
    }
    
    
    /**
     * Draw results in images
     */
    public void drawResults(ImagePlus img, Objects3DIntPopulation bactPop, String fileName, String imgName, String outDir) {
        ImageHandler imgObjects = ImageHandler.wrap(img).createSameDimensions();
        bactPop.drawInImage(imgObjects);
        IJ.run(imgObjects.getImagePlus(), "glasbey on dark", "");
        
        ImagePlus[] imgColors = {imgObjects.getImagePlus(), null, null, img};
        ImagePlus imgOut = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        imgOut.setCalibration(cal);
        FileSaver ImgObjectsFile = new FileSaver(imgOut);
        ImgObjectsFile.saveAsTiff(outDir+imgName+fileName);
        
        imgObjects.closeImagePlus();
        flush_close(imgOut);
    }
    
}
