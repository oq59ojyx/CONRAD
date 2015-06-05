package edu.stanford.rsl.tutorial.icecream;
import java.lang.Math;

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.*;
public class Phantom extends Grid2D {

	private float detectorSpacing;
	
	public Phantom(int size){
		super(size, size);
		createPhantom(size);
	}
	
	/*private void setDetectorSpacing(float ds) {
		this.detectorSpacing = ds;
	}
	
	public float getDetectorSpacing() {
		return detectorSpacing;
	} */
	
	private void createPhantom(int size){
		
		for (int i = 20;i < size; i++)
			for (int j = 20; j < size; j++){
				if ( (i*i+ j*j) < 12*size){
					this.setAtIndex(i, j, 15);
				}	
			}
		for (int i = size/2;i < size; i++)
			for (int j = size/2; j < size; j++){
				this.setAtIndex(i, j, 100);
			}
		for (int i = 0;i < 7*size/8; i++)
			this.setAtIndex(i, size/3, 50);

	}

	public Grid2D createSinogram(int numberProj, float detectorSpacing, int numberDetPixel, float SID){
		//setDetectorSpacing(detectorSpacing);
		Grid2D sinogram = new Grid2D(numberProj, numberDetPixel);
		for (int theta = 0; theta < numberProj; theta++)
		{
			double alpha =  ((2*Math.PI/(numberProj))*theta);

			for (int detPixel =0; detPixel< numberDetPixel; detPixel ++){

				double xspacePhantom = this.getSpacing()[0]/2;
			
				for (double t = -SID; t <  SID; t = t + xspacePhantom){
					//Umrechnung von Detektorkoord. in Weltkoord.
					float r = detPixel*detectorSpacing- (detectorSpacing*(numberDetPixel-1)/2);
					//Phantom
					float x = (float) ((float) r*Math.cos(alpha) + t* Math.sin(alpha)); 
					float y = (float) ((float) r*Math.sin(alpha) - t* Math.cos(alpha));

					 double[] id =  physicalToIndex(x, y);

					if (id[0] >= 0 && id[0] < this.getWidth()){
						if (id[1] >= 0 && id[1] < this.getHeight()){							
							sinogram.addAtIndex(theta, detPixel, InterpolationOperators.interpolateLinear(this, id[0], id[1]));
						}
					}
				}
			}
		}
		sinogram.show();
		return sinogram;
	}

	public Grid2D createFanogram(int numberProj, float detectorSpacing, int numberDetPixel, float dSD, int rotAngle, float dSI){
		
		Grid2D fanogram = new Grid2D(numberProj*numberDetPixel, numberProj);
		int idx = 0;
		for (int beta = 0; beta < numberProj; beta+= rotAngle)
		{
			//Umrechnung in RAD
			double alpha =  ((Math.PI*beta)/180);
			double sourcePos [] = {Math.sin(alpha)*dSI, Math.cos(alpha)*dSI};
			
			for (int detPixel =0; detPixel< numberDetPixel; detPixel ++)
			{
				float tMax = numberDetPixel*detectorSpacing/dSD;

				double xspacePhantom = this.getSpacing()[0]/2;
				double tdiff = tMax/ numberDetPixel;
				float t = (float) (tMax- (tdiff*(numberDetPixel-1)/2));
				t = t *detPixel -(detPixel/2-1);
				idx++;
				float ray = (float) (dSD/Math.cos(t));
				
					for (double l = 0; l< ray; l+= xspacePhantom)
					{
					
					//Umrechnung von Detektorkoord. in Weltkoord.
					
					//Phantom
					//float x =  (float) ((float) sourcePos[0]+ l* Math.cos(t * detPixel)); 
					//float y = (float) ((float)  sourcePos[1]+ l* Math.sin(t * detPixel));
						
						float gamma = (float) (-alpha - t);
						
						float x1 = (float) (l * Math.sin(gamma));
						float y1 = (float) (l * Math.cos(gamma));
						
						float x = (float) (x1- sourcePos[0]);
						float y = (float) (y1 - sourcePos[1]);
						
					 double[] id =  physicalToIndex(x, y);

					if (id[0] >= 0 && id[0] < this.getWidth()){
						if (id[1] >= 0 && id[1] < this.getHeight()){							
							fanogram.addAtIndex(idx, beta, InterpolationOperators.interpolateLinear(this, id[0], id[1]));
						}
					}
				}
			}
		}
		fanogram.show("fanogram");
		return fanogram;
	}
	
	public Grid2D rebinning(Grid2D fanogram, float dSI){
		
		Grid2D sinogram = new Grid2D (fanogram.getHeight(), fanogram.getWidth());
		
		for (int x = 0; x < sinogram.getHeight(); x++){
			for(int y = 0; y < sinogram.getWidth(); y ++){
				sinogram.addAtIndex(x, y, InterpolationOperators.interpolateLinear(fanogram, dSI * Math.sin(x), x+y));
			}
		}
		
		
		
		
		
		return sinogram;
	}
	

	public static void main(String[] args) {
		new ImageJ();
		int size = 256;
		Phantom p = new Phantom(size);
		p.setSpacing(0.1, 0.1);
		p.setOrigin(-(size-1)*p.spacing[0]/2, -(size-1)*p.spacing[1]/2 );
		p.show();	
		float min = NumericPointwiseOperators.min(p);
		float mean = NumericPointwiseOperators.mean(p);
		System.out.println(min);
		System.out.println(mean);
		
		float d = (float) (Math.sqrt(2)*p.getHeight()*p.getSpacing()[0]);
		float detectorSpacing = (float) 0.1;
		//p.createSinogram(180, detectorSpacing, (int) ((int) d/detectorSpacing), d/2 );
		p.createFanogram(20,  detectorSpacing, (int) ((int) d/detectorSpacing), d/2, 1, d/5);
		//p.createSinogram(180, (float)1.2, 500, d);
		
	}



}