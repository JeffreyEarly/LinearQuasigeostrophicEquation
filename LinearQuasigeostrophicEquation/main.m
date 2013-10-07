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
	    
		[GLVariable setPrefersSpatialMultiplication: NO];
		
		// Reasonable parameters to nondimensionalize by.
		GLFloat N_QG = 1.3; // cm
		GLFloat T_QG = 12; // days
		GLFloat L_QG = 47; // km
		
		NSLog(@"Define the dimensions of the problem.");
		
		// 256x128 at takes 11 second with optimized code.
		GLDimension *xDim = [[GLDimension alloc] initPeriodicDimension: YES nPoints: 256 domainMin: -1500/L_QG length: 2000/L_QG];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initPeriodicDimension: YES nPoints: 128 domainMin: -500/L_QG length: 1000/L_QG];
		yDim.name = @"y";
		GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: [NSArray arrayWithObject: [NSNumber numberWithDouble: 0.0]]];
		tDim.name = @"time";
		
		// Variables are always tied to a particular equation---so we create an equation object first.
		GLEquation *equation = [[GLEquation alloc] init];
		
		NSArray *spatialDimensions = [NSArray arrayWithObjects: xDim, yDim, nil];
		GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
		GLVariable *y = [GLVariable variableOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
		
		NSLog(@"Create and cache the differential operators that we will be using");
		
		// At the moment we know that this is the spectral operators, although in the future we'll have to set this up explicitly.
		GLSpectralDifferentialOperatorPool *diffOperators = [equation defaultDifferentialOperatorPoolForVariable: x];
		
		// Create the operator xx+yy-1---this is how you compute y from eta
		GLSpectralDifferentialOperator *laplacianMinusOne = [[diffOperators harmonicOperator] scalarAdd: -1.0];
		[diffOperators setDifferentialOperator: laplacianMinusOne forName: @"laplacianMinusOne"];
		
		// Create the operator 1/(xx+yy-1)---this is how you compute eta from y.
		GLSpectralDifferentialOperator *diffOp = [laplacianMinusOne scalarDivide: 1.0];
		[diffOperators setDifferentialOperator: diffOp forName: @"inverseLaplacianMinusOne"];
		
		// Because this equation is completely linear, we can make an operator that computes f from y directly
		GLSpectralDifferentialOperator *diffX = [diffOperators differentialOperatorWithName: @"x"];
		[diffOperators setDifferentialOperator: [[diffX negate] multiply: diffOp] forName: @"fFromY"];
		
		NSLog(@"Create the variable for the initial condition.");
		
		// Want, gaussian = amplitude * exp( - ((x-x0)*(x-x0) + (y-y0)*(y-y0))/(length*length) );
		GLFloat amplitude = 15.0/N_QG;
		GLFloat length = 80/L_QG;
		
		GLVariable *r2 = [[x times: x] plus: [y times: y]];
		GLVariable *gaussian = [[[r2 scalarMultiply: -1.0/(length*length)] exponentiate] scalarMultiply: amplitude];
		
		NSLog(@"Create a NetCDF to store the evolution of the gaussian");
		
		// Now we create a mutable variable in order to record the evolution of the Gaussian.
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: @"/Users/jearly/Desktop/LinearQuasigeostrophy.nc"] forEquation: equation overwriteExisting: YES];
		GLMutableVariable *sshHistory = [gaussian variableByAddingDimension: tDim];
		sshHistory.name = @"SSH";
		sshHistory = [netcdfFile addVariable: sshHistory];
		
		CGFloat cfl = 0.5;
		CGFloat deltaX = xDim.domainLength / ( (GLFloat) xDim.nPoints );
		GLFloat timeStep = cfl * deltaX / (L_QG/T_QG);
		GLFloat maxTime = 365/T_QG;
		
		NSLog(@"Time step the Gaussian with our equation: %d time steps", (int) (maxTime/timeStep));
		
		y = [gaussian diff: @"laplacianMinusOne"];
		GLIntegrationOperation *integrator = [GLIntegrationOperation rungeKutta4AdvanceY: y stepSize: timeStep fFromY:^(GLVariable *yNew) {
			return [yNew diff:@"fFromY"];
		}];
		
		for (GLFloat time = 0; time < maxTime; time += 1/T_QG)
		{
            @autoreleasepool {
				y = [integrator stepForward: y toTime: time];
				
				NSLog(@"Logging day: %f, step size: %f.", (integrator.currentTime*T_QG), integrator.lastStepSize*T_QG);
				// We're using spectral code, so it's possible (and is in fact the case) that the variable is not in the spatial domain.
				[tDim addPoint: [NSNumber numberWithDouble: integrator.currentTime]];
				GLVariable *eta = [[y diff: @"inverseLaplacianMinusOne"] spatialDomain];
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

