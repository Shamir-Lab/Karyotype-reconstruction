import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import ilog.concert.*;
import ilog.cplex.*;

public class ex1 {

	static int max_val = 20;	// Maximum value of an edge
	static String dir = ""; // This should be the path for the dir where the code is run from.
	
	// Params that are changed per simulation and should move out
	static String ver = ""; // This should be identical to the version string in the python code, 
							// it will be concatenated to all file names
	static double[] alpha_values = {0.1};	// LEGACY. For each sample the code iterates over all alpha values.

	
	public static void ILPRun(BPGraph g, String outfile, double alpha) throws FileNotFoundException, UnsupportedEncodingException
	{
		try {
			IloCplex cplex = new IloCplex(); 
			
			IloNumVar[][] xVarEdges = new IloIntVar[g.num_nodes][g.num_nodes];
			IloNumVar[][] xRefEdges = new IloIntVar[g.num_nodes][g.num_nodes];		
			IloNumVar[][] xIntEdges = new IloIntVar[g.num_nodes][g.num_nodes];		 
			
			// Create the vars - if an edge dosent exist the var is set to 0
			// The reason the ifs are not exclusive is because sometime ref edges and var edges can be parallel 
			for (int i=0; i<g.num_nodes; i++)
				for (int j=0; j<g.num_nodes; j++)
				{
					if (!g.IsEdge(i, j))
						xVarEdges[i][j] = cplex.intVar(0, 0);
					if (g.IsVarEdge(i, j))
						xVarEdges[i][j] = cplex.intVar(0, max_val);
					if (g.IsRefEdge(i, j))
						xRefEdges[i][j] = cplex.intVar(0, max_val);
					if (g.IsIntEdge(i, j))
						xIntEdges[i][j] = cplex.intVar(0, max_val);						
				}			

			// Add path constraints - the in and out quantity is equal
			for (int i=0; i<g.num_nodes; i++)
			{
				if (g.IsTerminusNode(i))
					continue;
				
				IloLinearNumExpr sumIntIn  = cplex.linearNumExpr();
				IloLinearNumExpr sumIntOut = cplex.linearNumExpr();
				IloLinearNumExpr sumRefIn  = cplex.linearNumExpr();
				IloLinearNumExpr sumRefOut = cplex.linearNumExpr();
				IloLinearNumExpr sumVarIn  = cplex.linearNumExpr();
				IloLinearNumExpr sumVarOut = cplex.linearNumExpr();
				
				for (int j=0; j<g.num_nodes; j++)
				{
					if (g.IsIntEdge(j, i))
						sumIntIn.addTerm(1, xIntEdges[j][i]);
					if (g.IsIntEdge(i, j))
						sumIntOut.addTerm(1, xIntEdges[i][j]);
					if (g.IsRefEdge(j, i))
						sumRefIn.addTerm(1, xRefEdges[j][i]);
					if (g.IsRefEdge(i, j))
						sumRefOut.addTerm(1, xRefEdges[i][j]);
					if (g.IsVarEdge(j, i))
						sumVarIn.addTerm(1, xVarEdges[j][i]);
					if (g.IsVarEdge(i, j))
						sumVarOut.addTerm(1, xVarEdges[i][j]);
				}
				
				cplex.addEq(sumIntIn, cplex.sum(sumRefOut, sumVarOut));
				cplex.addEq(sumIntOut, cplex.sum(sumRefIn, sumVarIn));
			}
						
			IloNumExpr[] yVarEdegs = new IloNumExpr[g.numUniqueVarEdges];	
			IloNumExpr[] yIntEdegs = new IloNumExpr[g.numIntEdegs/2];

			for (int i=0, var_idx = 0, int_idx = 0; i<g.num_nodes; i++)
			{
				for (int j=i; j<g.num_nodes; j++)
				{
					if (!g.IsEdge(i,j))
						continue;
					
					if (g.IsIntEdge(i, j))
					{
						yIntEdegs[int_idx] = cplex.prod(g.GetIntLen(i, j), cplex.abs(cplex.diff(cplex.sum(xIntEdges[i][j], xIntEdges[j][i]), g.GetIntEdgeWeight(i, j))));
						int_idx++;
					}
					if (g.IsVarEdge(i, j))
					{
						yVarEdegs[var_idx] = cplex.prod(g.GetVarEdgeWeight(i, j)*alpha, cplex.diff(1.0, cplex.min(1.0, cplex.sum(xVarEdges[i][j], xVarEdges[j][i]))));
						var_idx++;
					}
				}
			}
			cplex.addMinimize(cplex.sum(cplex.sum(yIntEdegs), cplex.sum(yVarEdegs)));
			
			
			if(cplex.solve())
			{
				PrintWriter writer = new PrintWriter(outfile, "UTF-8");
				PrintWriter sol_writer = new PrintWriter(outfile + "_sol_value.txt", "UTF-8");
				
				sol_writer.println("Solution value = " + cplex.getObjValue());
				
				writer.println(g.num_nodes);
				double int_edge_sum = 0;
				double var_edge_sum = 0;
				
				for (int i=0; i<g.num_nodes; i++)
				{
					for (int j=i; j<g.num_nodes; ++j)
					{
						if (g.IsVarEdge(i, j))
						{
							writer.println("v\t" + i + "\t" + j + "\t" + cplex.getValue(xVarEdges[i][j]) + "\t" + cplex.getValue(xVarEdges[j][i]) + "\t" + g.varEdgesWeights[i][j]);
							var_edge_sum += (1-Math.min(1, cplex.getValue(xVarEdges[i][j])+cplex.getValue(xVarEdges[j][i])))*g.varEdgesWeights[i][j]*alpha;
						}
						if (g.IsIntEdge(i,  j))
						{
							writer.println("i\t" + i + "\t" + j + "\t" + cplex.getValue(xIntEdges[i][j]) + "\t" + cplex.getValue(xIntEdges[j][i]) + "\t" + g.intEdgesWeights[i][j]);
							int_edge_sum += Math.abs((cplex.getValue(xIntEdges[i][j])+cplex.getValue(xIntEdges[j][i]))-g.intEdgesWeights[i][j]) * g.intEdgesLengths[i][j];
						}
						if (g.IsRefEdge(i, j))
						{
							writer.println("r\t" + i + "\t" + j + "\t" + cplex.getValue(xRefEdges[i][j]) + "\t" + cplex.getValue(xRefEdges[j][i]));
						}
					}
				}
				sol_writer.println("int sum = " + int_edge_sum);
				sol_writer.println("var sum = " + var_edge_sum);
				
				writer.close();
				sol_writer.close();
			}
			
			cplex.end();
		}
		
		catch (IloException e) {
			System.err.println("Concert exception '" + e + "' caught");
		}
	}
	
	public static void main(String[] args) throws IOException {

		String samples_file = args[0];
		BPGraph g;
		String filename;
		

		 //Create object of FileReader
	    FileReader inputFile = new FileReader(dir+samples_file);

	    //Instantiate the BufferedReader Class
	    BufferedReader bufferReader = new BufferedReader(inputFile);

	    String line;	    
	    String delims = "[\t]";
	    String sample_name;
	    // Read file line by line and print on the console
	    while ((line = bufferReader.readLine()) != null)
	    {
	    	// Here build the graph
	    	String[] tokens = line.split(delims);
	    	sample_name = tokens[0] +"_" + tokens[1];
	    	
	    	filename = dir + "/Results/" + tokens[0] + "/" + sample_name + "/" + sample_name;
	    	for (int i=0; i<alpha_values.length; i++)
			{
				g = new BPGraph(filename + "_graph.txt");
				
				System.out.println("Processing sample " + sample_name + "_ILP_v"+ver +"_" + alpha_values[i]);
				ILPRun(g, filename + "_ILP_v" + ver +"_alpha" + alpha_values[i] + ".txt", alpha_values[i]);
			}
	    }
	    System.out.println("Done.");
	}
}
