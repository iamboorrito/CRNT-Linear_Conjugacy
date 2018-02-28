/*********************************************
 * OPL 12.8.0.0 Model
 * Author: Evan Burton
 * Creation Date: Feb 25, 2018 at 6:31:22 PM
 *********************************************/

// Clocking
/////////////////////////
float temp;
execute{
var before = new Date();
temp = before.getTime();
}
/////////////////////////

// Data declarations
int rows = ...;
int cols = ...;
float Y[1..rows][1..cols] = ...;
float M[1..rows][1..cols] = ...;
int ubound = ...;
float eps = ...;

// Coordinate pair data structure needed for iteration over
// custom sets of indices
tuple pair {
	int i;
	int j;
};

// This is the set of pairs (i,j) s.t. that i != j
{pair} off_diag = {<i, j> | i,j in 1..cols: i != j};

// Decision variable declarations
dvar float Ab[1..cols][1..cols];
dvar float T[1..rows];
dvar float Ah[1..cols][1..cols];
// Should come up with better way of storing/indexing delta since there are
// only n(n-1) values used and n ignored.
dvar int delta[1..cols][1..cols];
dvar int z;

// Objective function
minimize z;

// Begin constraint declarations
subject to {

	z == sum(<i,j> in off_diag) delta[i][j];
	z == 13;

	// Y Ab = T M
	forall(i in 1..rows, j in 1..cols){
		sum(k in 1..cols) Y[i][k]*Ab[k][j] - T[i]*M[i][j] == 0;
	}
	
	// Constraints on sums of entries in Ab, Ah
	forall(j in 1..cols){
	
		// Columns of Ab sum to 0	
		sum(i in 1..cols) Ab[i][j] == 0;
		// Ah weak reversibility constraints
		sum(i in 1..cols) Ah[i][j] == 0;
		sum(i in 1..cols) Ah[j][i] == 0;
	}
	
	// Constraints on the off diagonal entries of Ab and Ah
	forall(<i,j> in off_diag){
			// Ab constraints
			Ab[i][j] >= 0;	
			Ab[i][j] - eps*delta[i][j] >= 0;
		    Ab[i][j] - ubound*delta[i][j] <= 0;
		   
		   	// Ah constraints
		   	Ah[i][j] >= 0;	
			Ah[i][j] - eps*delta[i][j] >= 0;
		    Ah[i][j] - ubound*delta[i][j] <= 0;
		   
		    // delta is binary
		    0 <= delta[i][j] <= 1;
	}
	
	// Bound diagonal entries of Ab and Ah
	forall(j in 1..cols){	
	
			Ab[j][j] <= 0;
			
			Ah[j][j] <= 0;	
			
	}
	
	// Bounds on T
	forall(i in 1..rows){
		eps <= T[i] <= 1/eps;	
	}

}

// Clocking
//////////////////////////
execute{
var after = new Date();
writeln("solving time ~= ", (after.getTime()-temp));
}
//////////////////////////

execute DISPLAY{};

	
	/*
	// Old code
	// Bounds on entries of Ab, Ak, and delta binary
	forall(i in 1..cols, j in 1..cols){
		if(i != j){
			// Ab constraints
			Ab[i][j] >= 0;	
			Ab[i][j] - eps*delta[i][j] >= 0;
		    Ab[i][j] - ubound*delta[i][j] <= 0;
		   
		   	// Ah constraints
		   	Ah[i][j] >= 0;	
			Ah[i][j] - eps*delta[i][j] >= 0;
		    Ah[i][j] - ubound*delta[i][j] <= 0;
		   
		    // delta is binary
		    0 <= delta[i][j] <= 1;
		}else{
			Ab[i][j] <= 0;
			Ab[i][j] + ubound*delta[i][j] >= 0;
			
			Ah[i][j] <= 0;	
			Ah[i][j] + ubound*delta[i][j] >= 0;
		}	
	}
	*/