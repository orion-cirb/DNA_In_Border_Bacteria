package FociBacteria_Tools;

import FociBacteria_Tools.Cellpose.CellposeTaskSettings;
import FociBacteria_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import fiji.util.gui.GenericDialogPlus;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.process.AutoThresholder;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Voxel3D;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.VoxelInt;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.geom2.measurements.MeasureFeret;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationClosestDistance;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationDistance;
import mcib3d.geom2.measurementsPopulation.PairObjects3DInt;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import org.apache.commons.io.FilenameUtils;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;


/**
 * @author Orion-CIRB
 */
public class Tools {
    private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
      
    public Calibration cal = new Calibration();
    private double pixelSurf = 0;
    String[] channelsName = {"Bacteria: ", "Foci ch1: ", "Foci ch2: "};
    
     // Omnipose
    private String omniposeEnvDirPath = "/opt/miniconda3/envs/omnipose";
    private String omniposeModelsPath = System.getProperty("user.home")+"/.cellpose/models/";
    private String omniposeModel = "bact_phase_omnitorch_0";
    private int omniposeDiameter = 10;
    private int omniposeMaskThreshold = 0;
    private double omniposeFlowThreshold = 0.4;
    private boolean useGpu = true;
    
    // Bacteria
    private double minBactSurface = 1;
    private double maxBactSurface = 10;
    
    // Foci segmentation method
    private double fociDOGMin = 1;
    private double fociDOGMax = 2;
    private String[] fociThMethods = AutoThresholder.getMethods();
    public String foci1Th = fociThMethods[6];
    public String foci2Th = fociThMethods[11];
    
    // Foci
    private double minFociSurface = 0.01;
    private double maxFociSurface = 0.5;

    private final CLIJ2 clij2 = CLIJ2.getInstance();
    
    
    
    
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
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.showMessage("Error", "CLIJ not installed, please install from update site");
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
        for (String ch : channels) {
            gd.addChoice(channelsName[index], channels, ch);
            index++;
        }
        
        gd.addMessage("Bacteria detection", Font.getFont("Monospace"), Color.blue);
        if (IJ.isWindows()) {
            omniposeEnvDirPath = System.getProperty("user.home")+"\\miniconda3\\envs\\omnipose";
            omniposeModelsPath = System.getProperty("user.home")+"\\.cellpose\\models\\";
        }
        gd.addDirectoryField("Omnipose environment directory: ", omniposeEnvDirPath);
        gd.addDirectoryField("Omnipose models path: ", omniposeModelsPath);
        gd.addNumericField("Min bacterium surface (µm2): ", minBactSurface);
        gd.addNumericField("Max bacterium surface (µm2): ", maxBactSurface);
        
        gd.addMessage("Foci detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min foci surface (µm2): ", minFociSurface);
        gd.addNumericField("Max foci surface (µm2): ", maxFociSurface);
        gd.addChoice("Foci (CFP)    :",fociThMethods, foci1Th);
        gd.addChoice("Foci (phiYFP) :",fociThMethods, foci2Th);
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm):", cal.pixelWidth);
        gd.showDialog();
        
        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();
        if(gd.wasCanceled())
           ch = null;
        
        omniposeEnvDirPath = gd.getNextString();
        omniposeModelsPath = gd.getNextString();
        minBactSurface = (float) gd.getNextNumber();
        maxBactSurface = (float) gd.getNextNumber();
        
        minFociSurface = (float) gd.getNextNumber();
        maxFociSurface = (float) gd.getNextNumber();
        foci1Th = gd.getNextChoice();
        foci2Th = gd.getNextChoice();
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = 1;
        pixelSurf = cal.pixelWidth*cal.pixelHeight;
        
        return(ch);
    }
    
    
    /**
     * Do Z projection
     */
    public ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    
    /**
    * Detect bacteria with Omnipose
    */
    public Objects3DIntPopulation omniposeDetection(ImagePlus imgBact){
        ImagePlus imgIn = new Duplicator().run(imgBact);
        imgIn.setCalibration(cal);
        
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
        PrintStream console = System.out;
        System.setOut(new NullPrintStream());
        ImagePlus imgOut = cellpose.run();
        System.setOut(console);
        imgOut.setCalibration(cal);
        
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
     * Find foci with DOG
     */
    public Objects3DIntPopulation findFoci(ImagePlus img, String fociTh) {
        ImagePlus imgFilter = DOG(img, fociDOGMin, fociDOGMax);
        ImagePlus imgBin = threshold(imgFilter, fociTh);
        imgBin.setCalibration(cal);
        
        Objects3DIntPopulation fociPop = getPopFromImage(imgBin);
        fociPop = new Objects3DIntPopulationComputation(fociPop).getFilterSize(minFociSurface/pixelSurf, maxFociSurface/pixelSurf);
        
        flush_close(imgFilter);
        flush_close(imgBin);
        
        return(fociPop);
    }
    
    
    /**
     * Difference of gaussian
     */ 
    public ImagePlus DOG(ImagePlus img, double size1, double size2) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian2D(imgCL, imgCLDOG, size1, size1, size2, size2);
        ImagePlus imgDOG = clij2.pull(imgCLDOG);
        clij2.release(imgCL);
        clij2.release(imgCLDOG);
        return(imgDOG);
    }
    
    
    /**
     * Laplacian of gaussians
     */ 
    /*public ImagePlus LOG(ImagePlus img, double size) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer temp1 = clij2.create(imgCL);
        ClearCLBuffer imgCLLOG = clij2.create(imgCL);
        clij2.gaussianBlur2D(imgCL, temp1, size, size);
        clij2.laplaceBox(temp1, imgCLLOG);
        clij2.release(imgCL);
        clij2.release(temp1);
        ImagePlus imgLog = clij2.pull(imgCLLOG);
        clij2.release(imgCLLOG);
        return(imgLog);
    }*/
    
    
    /**
     * Threshold
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        clij2.release(imgCL);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
    /**
     * Get population form labelled image
     */
    public Objects3DIntPopulation getPopFromImage(ImagePlus img) {
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        return pop;
    }
    
    
     /**
     * Find foci in bacteria
     * Set label of mother bacterium in foci object
     */
    public void fociBactLink(Objects3DIntPopulation bactPop, Objects3DIntPopulation fociPop) {
        if (bactPop.getNbObjects() != 0 && fociPop.getNbObjects() != 0) {
            for (Object3DInt bact : bactPop.getObjects3DInt()) {
                for (Object3DInt foci : fociPop.getObjects3DInt()) {
                    MeasureCentroid fociCenter = new MeasureCentroid(foci);
                    if (bact.contains(fociCenter.getCentroidRoundedAsVoxelInt())){
                       foci.setIdObject(bact.getLabel()); 
                    }
                }
            }
        }
        // Remove foci not in bacteria
        fociPop.getObjects3DInt().removeIf(p -> p.getIdObject() == 0);
        fociPop.resetLabels();
    }
    
    
    /**
     * Compute bacteria parameters and save them in file
     * @throws java.io.IOException
     */
    public void saveResults(Objects3DIntPopulation bactPop, Objects3DIntPopulation foci1Pop, Objects3DIntPopulation foci2Pop, String imgName, BufferedWriter distFile, BufferedWriter colocFile) throws IOException {
        for (Object3DInt bact : bactPop.getObjects3DInt()) {
            float bactLabel = bact.getLabel();
            double bactSurf = new MeasureVolume(bact).getValueMeasurement(MeasureVolume.VOLUME_UNIT);
            VoxelInt feret1Unit = new MeasureFeret(bact).getFeret1Unit();
            VoxelInt feret2Unit = new MeasureFeret(bact).getFeret2Unit();
            double bactLength = feret1Unit.distance(feret2Unit)*cal.pixelWidth;
            
            Objects3DIntPopulation foci1BactPop = findFociInBact(bactLabel, foci1Pop);
            Objects3DIntPopulation foci2BactPop = findFociInBact(bactLabel, foci2Pop);
            MeasurePopulationDistance distances = new MeasurePopulationDistance(foci1BactPop, foci2BactPop);
            for(Object3DInt foci1: foci1BactPop.getObjects3DInt()) {
                Voxel3D foci1Center = new MeasureCentroid(foci1).getCentroidAsVoxel();
                double poleDist = feret1Unit.distance(foci1Center)*cal.pixelWidth;
                
                distFile.write(imgName+"\t"+bactLabel+"\t"+bactSurf+"\t"+bactLength+"\t"+foci1BactPop.getNbObjects()+"\t"+
                        foci2BactPop.getNbObjects()+"\t1\t"+foci1.getLabel()+"\t"+poleDist);
                for(Object3DInt foci2: foci2BactPop.getObjects3DInt()) {
                    distFile.write("\t" + distances.getValueObjectsPair(foci1, foci2)*cal.pixelWidth);
                }
                distFile.write("\n");
            }
            
            for(Object3DInt foci2: foci2BactPop.getObjects3DInt()) {
                Voxel3D foci2Center = new MeasureCentroid(foci2).getCentroidAsVoxel();
                double poleDist = feret1Unit.distance(foci2Center)*cal.pixelWidth;
                
                distFile.write(imgName+"\t"+bactLabel+"\t"+bactSurf+"\t"+bactLength+"\t"+foci1BactPop.getNbObjects()+"\t"+
                        foci2BactPop.getNbObjects()+"\t2\t"+foci2.getLabel()+"\t"+poleDist);
                for(Object3DInt foci1: foci1BactPop.getObjects3DInt()) {
                    distFile.write("\t" + distances.getValueObjectsPair(foci1, foci2)*cal.pixelWidth);
                }
                distFile.write("\n");
            }
            distFile.flush();
            
            boolean colocEvent = false;
            MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(foci1BactPop, foci2BactPop);
            for (Object3DInt foci1: foci1BactPop.getObjects3DInt()) {
                for (Object3DInt foci2: foci2BactPop.getObjects3DInt()) {
                    double colocVal = coloc.getValueObjectsPair(foci1, foci2);
                    if (colocVal > 0) {
                        colocEvent = true;
                        
                        Voxel3D foci1Center = new MeasureCentroid(foci1).getCentroidAsVoxel();
                        double poleDist1 = feret1Unit.distance(foci1Center)*cal.pixelWidth;
                        Voxel3D foci2Center = new MeasureCentroid(foci2).getCentroidAsVoxel();
                        double poleDist2 = feret1Unit.distance(foci2Center)*cal.pixelWidth;
                        
                        colocFile.write(imgName+"\t"+bactLabel+"\t"+bactSurf+"\t"+bactLength+"\t"+
                                foci1BactPop.getNbObjects()+"\t"+foci2BactPop.getNbObjects()+"\tYes\t"+
                                foci1.getLabel()+"\t"+poleDist1+"\t"+foci2.getLabel()+"\t"+poleDist2+"\n");
                    }
                }
            }
            if (!colocEvent) {
                colocFile.write(imgName+"\t"+bactLabel+"\t"+bactSurf+"\t"+bactLength+"\t"+
                                foci1BactPop.getNbObjects()+"\t"+foci2BactPop.getNbObjects()+"\tNo\n");
            }
            colocFile.flush();
            

            /*
            header = "Image name\t# bacterium\tBacterium surface (µm2)\tBacterium length (µm)\t" +
                     "Foci1 number\tFoci2 number\tFoci1 min. distance to bacterium\t" +
                     "Foci2 min. distance to bacterium\tFoci1-foci2 min distance\tFoci1-foci2 max distance\n"
            
            Objects3DIntPopulation foci1BactPop = findFociBact(bactLabel, foci1Pop);
            int foci1Nb = foci1BactPop.getNbObjects();
            ArrayList<Double> minFoci1Dist = (foci1Nb == 0) ? null : fociBactDistance(foci1BactPop, feret1Unit);
            Objects3DIntPopulation foci2BactPop = findFociBact(bactLabel, foci2Pop);  
            int foci2Nb = foci2BactPop.getNbObjects();
            ArrayList<Double> minFoci2Dist = (foci2Nb == 0) ? null : fociBactDistance(foci2BactPop, feret1Unit);
            int fociMax = Math.max(foci1Nb, foci2Nb);
            double[] foci1Foci2Distance = (foci1Nb == 0 || foci2Nb == 0) ? null : fociFociDistance(foci1BactPop, foci2BactPop);
            file.write(imgName+"\t"+bactLabel+"\t"+bactSurf+"\t"+bactLength+"\t"+foci1Nb+"\t"+foci2Nb+"\t");
            if (fociMax != 0) {
                for (int i = 0; i < fociMax; i++) {
                    if (i != 0)
                        file.write("\t\t\t\t\t\t");
                    if (minFoci1Dist != null && i < minFoci1Dist.size())
                        file.write(minFoci1Dist.get(i)+"\t");
                    else
                        file.write("\t");
                    if (minFoci2Dist != null && i < minFoci2Dist.size())
                        file.write(minFoci2Dist.get(i)+"\t");
                    else
                        file.write("\t");
                    if (foci1Foci2Distance != null && i == 0)
                        file.write(foci1Foci2Distance[0]+"\t"+foci1Foci2Distance[1]+"\n");
                    else
                        file.write("\t\t\n");
                }
            }
            else
                file.write("\t\t\t\t\n");*/
        }
    }
    
    
    /**
     * Get foci in bacterium
     */
    private Objects3DIntPopulation findFociInBact(float bactLabel, Objects3DIntPopulation fociPop) {
        Objects3DIntPopulation fociBactPop = new Objects3DIntPopulation();
        for (Object3DInt foci: fociPop.getObjects3DInt()) {
                if (foci.getIdObject() == bactLabel)
                    fociBactPop.addObject(foci);
        }
        fociBactPop.resetLabels();
        return(fociBactPop);           
    }
    
    
    /**
     * Compute bacterium area
     */
    private double bacteriumSurface(Object3DInt bactObj) {
        int pixelNb = bactObj.getObject3DPlanes().get(0).getVoxels().size();
        return(pixelNb*pixelSurf);
    }
    
    
    /**
     * Compute min distance between foci to bacteria border
     */
    private ArrayList<Double> fociBactDistance(Objects3DIntPopulation fociPop, VoxelInt bactFeret1) {
        ArrayList<Double> minFociDist = new ArrayList<>();
        for (Object3DInt foci : fociPop.getObjects3DInt()) {
            Voxel3D fociCenter = new MeasureCentroid(foci).getCentroidAsVoxel();
            double dist = bactFeret1.distance(fociCenter)*cal.pixelWidth;
            minFociDist.add(dist);
        }
        return(minFociDist);
    }
    
    
    /**
     * Compute min,max distances between foci1 and foci2
     */
    private double[] fociFociDistance(Objects3DIntPopulation foci1Pop, Objects3DIntPopulation foci2Pop) {
        MeasurePopulationClosestDistance allDist =  new MeasurePopulationClosestDistance​(foci1Pop, foci2Pop);
        allDist.setDistanceMax(5);
        List<PairObjects3DInt> allDistPairs = allDist.getAllPairs(true);
        double[] minMaxDist = {allDistPairs.get(0).getPairValue(), allDistPairs.get(allDistPairs.size()-1).getPairValue()};
        return(minMaxDist);
    }
    
    
    // Save objects image
    public void drawResults(ImagePlus img, Objects3DIntPopulation bactPop, Objects3DIntPopulation foci1Pop, Objects3DIntPopulation foci2Pop, String imgName, String outDir) {
        ImageHandler imgBact = ImageHandler.wrap(img).createSameDimensions();
        bactPop.drawInImage(imgBact);
        IJ.run(imgBact.getImagePlus(), "glasbey on dark", "");
        ImagePlus[] imgColors1 = {imgBact.getImagePlus(), null, null, img};
        ImagePlus imgOut1 = new RGBStackMerge().mergeHyperstacks(imgColors1, true);
        imgOut1.setCalibration(cal);
        FileSaver ImgObjectsFile1 = new FileSaver(imgOut1);
        ImgObjectsFile1.saveAsTiff(outDir+imgName+"_bacteria.tif");
        
        ImageHandler imgFoci1 = ImageHandler.wrap(img).createSameDimensions();
        foci1Pop.drawInImage(imgFoci1);
        ImageHandler imgFoci2 = ImageHandler.wrap(img).createSameDimensions();
        foci2Pop.drawInImage(imgFoci2);
        ImagePlus[] imgColors2 = {imgFoci1.getImagePlus(), imgFoci2.getImagePlus(), null, img};
        ImagePlus imgOut2 = new RGBStackMerge().mergeHyperstacks(imgColors2, false);
        imgOut2.setCalibration(cal);
        FileSaver ImgObjectsFile2 = new FileSaver(imgOut2);
        ImgObjectsFile2.saveAsTiff(outDir + imgName + "_foci.tif");
        
        flush_close(imgBact.getImagePlus());
        flush_close(imgFoci1.getImagePlus());
        flush_close(imgFoci2.getImagePlus());
        flush_close(imgOut1);
        flush_close(imgOut2);
    }
    
}
