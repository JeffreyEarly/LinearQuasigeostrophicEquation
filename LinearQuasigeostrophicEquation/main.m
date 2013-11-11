//
//  main.m
//  LinearQuasigeostrophicEquation
//
//  Created by Jeffrey Early on 3/21/12.
//  Copyright (c) 2012 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

int main (int argc, const char * argv[])
{
	
	@autoreleasepool {
		
		// Reasonable parameters to nondimensionalize by.
		GLFloat N_QG = 1.3; // cm
		GLFloat T_QG = 12; // days
		GLFloat L_QG = 47; // km
		
		/************************************************************************************************/
		/*		Define the problem dimensions															*/
		/************************************************************************************************/
		
		// 256x128 at takes 11 second with optimized code.
		GLDimension *xDim = [[GLDimension alloc] initPeriodicDimension: YES nPoints: 256 domainMin: -1500/L_QG length: 2000/L_QG];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initPeriodicDimension: YES nPoints: 128 domainMin: -500/L_QG length: 1000/L_QG];
		yDim.name = @"y";
		GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
		tDim.name = @"time";
		
		// Variables are always tied to a particular equation---so we create an equation object first.
		GLEquation *equation = [[GLEquation alloc] init];
		
		NSArray *spatialDimensions = @[xDim, yDim];
		GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
		GLVariable *y = [GLVariable variableOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
		
		/************************************************************************************************/
		/*		Create and cache the differential operators we will need								*/
		/************************************************************************************************/
		
		NSArray *spectralDimensions = [x dimensionsTransformedToBasis: x.differentiationBasis];
		
		GLLinearTransform *laplacian = [GLLinearTransform harmonicOperatorFromDimensions: spectralDimensions forEquation: equation];
		GLLinearTransform *laplacianMinusOne = [laplacian plus: @(-1.0)];
		GLLinearTransform *inverseLaplacianMinusOne = [laplacianMinusOne inverse];
        
		GLLinearTransform *diffX = [GLLinearTransform differentialOperatorWithDerivatives:@[@(1),@(0)] fromDimensions:spectralDimensions forEquation:equation];
		GLLinearTransform *fFromY = [[diffX times: @(-1)] times: inverseLaplacianMinusOne];
		
        
        //GLLinearTransform *diffOp = [GLLinearTransform differentialOperatorWithDerivatives: @[@(0),@(1)] fromDimensions: spectralDimensions forEquation: equation];
        //GLLinearTransform *diffOp = [GLLinearTransform differentialOperatorOfOrder: 1 fromDimension: spectralDimensions.firstObject forEquation: equation];
//        GLLinearTransform *diffOp = fFromY;
//        
//        [diffOp solve];
//        GLFloat *val = diffOp.pointerValue;
//        NSUInteger n = [diffOp.fromDimensions.lastObject nPoints];
//        for (NSUInteger i=0; i<diffOp.nDataPoints; i++) {
//            if (i%n==0) printf("\n");
//            printf("%g ", val[i]);
//        }
//		
//        return 0;
        
		/************************************************************************************************/
		/*		Create the initial conditions															*/
		/************************************************************************************************/
		
		// Want, gaussian = amplitude * exp( - ((x-x0)*(x-x0) + (y-y0)*(y-y0))/(length*length) );
		GLFloat amplitude = 15.0/N_QG;
		GLFloat length = 80/L_QG;
		
		GLVariable *r2 = [[x times: x] plus: [y times: y]];
		GLVariable *gaussian = [[[r2 times: @(-1.0/(length*length))] exponentiate] times: @(amplitude)];
        
		/************************************************************************************************/
		/*		Create a file to output data															*/
		/************************************************************************************************/
		
		// Now we create a mutable variable in order to record the evolution of the Gaussian.
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: @"/Users/jearly/Desktop/LinearQuasigeostrophy.nc"] forEquation: equation overwriteExisting: YES];
		GLMutableVariable *sshHistory = [gaussian variableByAddingDimension: tDim];
		sshHistory.name = @"SSH";
		sshHistory = [netcdfFile addVariable: sshHistory];
		
		/************************************************************************************************/
		/*		Estimate the time step size																*/
		/************************************************************************************************/
		
		CGFloat cfl = 0.5;
		CGFloat deltaX = xDim.domainLength / ( (GLFloat) xDim.nPoints );
		GLFloat timeStep = cfl * deltaX / (L_QG/T_QG);
		GLFloat maxTime = 365/T_QG;
		
		NSLog(@"Time step the Gaussian with our equation: %d time steps", (int) (maxTime/timeStep));
		
		/************************************************************************************************/
		/*		Create an integrator: dy/dt=f															*/
		/************************************************************************************************/
		
        y = [gaussian differentiateWithOperator: laplacianMinusOne];
		GLRungeKuttaOperation *integrator = [GLRungeKuttaOperation rungeKutta4AdvanceY: @[y] stepSize: timeStep fFromTY:^(GLScalar *t, NSArray *yNew) {
			return @[[fFromY transform: yNew[0]]];
		}];
		
		/************************************************************************************************/
		/*		Now iterate! Stop every day to write out some data.										*/
		/************************************************************************************************/
		
		for (GLFloat time = 0; time < maxTime; time += 1/T_QG)
		{
            @autoreleasepool {
				NSArray *yout = [integrator stepForwardToTime: time];
				
				NSLog(@"Logging day: %f, step size: %f.", (integrator.currentTime*T_QG), integrator.lastStepSize*T_QG);
				// We're using spectral code, so it's possible (and is in fact the case) that the variable is not in the spatial domain.
				[tDim addPoint: @(integrator.currentTime)];
				GLVariable *eta = [[inverseLaplacianMinusOne transform: yout[0]] spatialDomain];
                
                // eta is actually fine here, it's just that writing to the netcdf file is somehow broken.
				[sshHistory concatenateWithLowerDimensionalVariable: eta alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
            }
		}
		
		NSLog(@"Close the NetCDF file and wrap up");
		
		[equation waitUntilAllOperationsAreFinished];
		
		// The NetCDF file may still be writing data. We need to make sure it finishes before we exit.
		[netcdfFile waitUntilAllOperationsAreFinished];
		[netcdfFile close];
		
	    
	}
    return 0;
}

