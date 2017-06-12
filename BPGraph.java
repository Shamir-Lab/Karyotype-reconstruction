import java.io.*;

/**
 * This holds a bridge graph. A bridge graph contains three boolean matrices of size N X N representing
 * the three types of edges.
 * It also holds the number of nodes (N) and a counter for all types of edges.
 */
public class BPGraph {
	int num_nodes;
	int numVarEdegs;
	int numIntEdegs;
	int numRefEdegs;
	
	
	boolean[][] intEdges;
	boolean[][] refEdges;
	boolean[][] varEdges;
	
	// This is used to hold all the terminal edges in the graph and also to force the solution to be in one direction
	boolean[][] terminalEdges;
	
	double[][] intEdgesWeights;
	double[][] intEdgesLengths;		// This is holding the normalized lengths of the intervals
	double[][] refEdgesWeights;
	double[][] varEdgesWeights;		// Holding the normalized support according to the total support value
	
	
	public BPGraph(int n)
	{
		num_nodes = n;
		intEdges = new boolean[n][n];
		refEdges = new boolean[n][n];
		varEdges = new boolean[n][n];
		
		intEdgesWeights = new double[n][n];
		refEdgesWeights = new double[n][n];		// Not used for now
		varEdgesWeights = new double[n][n];
		intEdgesLengths = new double[n][n];
		
		terminalEdges = new boolean[n][n];
		
		numIntEdegs = 0;
		numRefEdegs = 0;
		numVarEdegs = 0;
	}
	
	public double GetVarEdgeWeight(int u, int v)
	{
		if (!IsVarEdge(u,v))
			return -1;
		return varEdgesWeights[u][v];
	}
	
	public double GetIntEdgeWeight(int u, int v)
	{
		if (!IsIntEdge(u,v))
			return -1;
		return intEdgesWeights[u][v];
	}
	
	public double GetIntLen(int u, int v)
	{
		if (!IsIntEdge(u,v))
			return -1;
		return intEdgesLengths[u][v];
	}
	
	public double GetRefEdgeWeight(int u, int v)
	{
		if (!IsRefEdge(u,v))
			return -1;
		return refEdgesWeights[u][v];
	}
	 
	public void SetRefEdge(int u, int v, double w)
	{		
		if (!refEdges[u][v])
			numRefEdegs++;
	
		refEdges[u][v] = true;
		refEdgesWeights[u][v] = w;
	}
	
	public void SetVarEdge(int u, int v, double w)
	{	
		if (u==v)
			return;
		if (!varEdges[u][v])
			numVarEdegs++;
	
		varEdges[u][v] = true;
		varEdgesWeights[u][v] = w;
	}
	
	public void SetIntEdge(int u, int v, double w, double length, boolean force_dir)
	{		
		if (!intEdges[u][v])
			numIntEdegs++;
	
		intEdges[u][v] = true;
		intEdgesWeights[u][v] = w;
		intEdgesLengths[u][v] = length;
		
		terminalEdges[u][v] = force_dir;
	}
	
	public int NumEdges()
	{
		return numRefEdegs + numIntEdegs + numVarEdegs;
	}
	
	public boolean IsEdge(int u, int v)
	{
		return (u<num_nodes && (intEdges[u][v] || refEdges[u][v] || varEdges[u][v]));
	}
	
	public boolean IsVarEdge(int u, int v)
	{
		return (u<num_nodes &&  varEdges[u][v]);
	}
	
	public boolean IsIntEdge(int u, int v)
	{
		return (u<num_nodes &&  intEdges[u][v]);
	}
	
	public boolean IsRefEdge(int u, int v)
	{
		return (u<num_nodes && refEdges[u][v]);
	}
	
	/** 
	 * Return true if this node is at the edge of a chromosme
	 */
	public boolean IsTerminusNode(int u)
	{
		int edge_counter=0;
		for (int i=0; i<num_nodes; i++)
		{
			if (IsRefEdge(i, u))
				return false;
			if (IsVarEdge(i, u))
				return false;
			if (IsRefEdge(u, i))
				return false;
			if (IsVarEdge(u, i))
				return false;
			// Really just a check for sanity of the graph
			if (IsIntEdge(i, u) || IsIntEdge(u, i))
				edge_counter++;
		}
		if (edge_counter!=1)
			System.out.println("graph has errors. node " + u + " has " + edge_counter + "int edges");
		return true;
	}
	
	/**
	 * Constructor fro, a text file.
	 * 
	 * @param fileName - Name of a file holding the information of the graph.
	 */
	public BPGraph(String fileName)
	{
		String delims = "[\t]";
		
		char e;
    	int u;
    	int v;
		
		System.out.println("Reading File from Java code");
	    
	    try{

	    //Create object of FileReader
	    FileReader inputFile = new FileReader(fileName);

	    //Instantiate the BufferedReader Class
	    BufferedReader bufferReader = new BufferedReader(inputFile);

	    //Variable to hold the one line data
	    String line = bufferReader.readLine();
	    
	    int graph_size = Integer.parseInt(line);
	    if (graph_size <= 0)	//Add more sanity checks
	    	return;		// Add a println
    
    	num_nodes = graph_size;
		intEdges = new boolean[graph_size][graph_size];
		refEdges = new boolean[graph_size][graph_size];
		varEdges = new boolean[graph_size][graph_size];
		
		intEdgesWeights = new double[graph_size][graph_size];
		refEdgesWeights = new double[graph_size][graph_size];
		varEdgesWeights = new double[graph_size][graph_size];
		intEdgesLengths = new double[graph_size][graph_size];
		
		terminalEdges = new boolean[graph_size][graph_size];

	    // Read file line by line and print on the console
	    while ((line = bufferReader.readLine()) != null)
	    {
	    	// Here build the graph
	    	String[] tokens = line.split(delims);
	    	e = tokens[0].charAt(0);
	    	u = Integer.parseInt(tokens[1]);
	    	v = Integer.parseInt(tokens[2]);
	    	
	    	// Each edge is represented as two directed edges
	    	if (e=='i')
	    	{
	    		double weight = Double.parseDouble(tokens[3]);
	    		double len = Double.parseDouble(tokens[4]);
	    		
	    		SetIntEdge(u, v, weight, len, false);
	    		SetIntEdge(v, u, weight, len, false);
	    		
	    		
	    		if (tokens.length > 5)
	    		{
	    			if (tokens[5].equals("->"))
	    				SetIntEdge(u, v, weight, len, true);
	    			else if (tokens[5].equals("<-"))
	    				SetIntEdge(v, u, weight, len, true);
	    		}
	    	}
	    	else if (e=='v')
	    	{
	    		double supp = Double.parseDouble(tokens[3]);
	    		SetVarEdge(u, v, supp);
	    		SetVarEdge(v, u, supp);
	    	}
	    	else if (e=='r')
	    	{
	    		SetRefEdge(u, v, 0);
	    		SetRefEdge(v, u, 0);
	    	}
	    }
	    //Close the buffer reader
	    bufferReader.close();
	    
	    }catch(Exception ex){
	            System.out.println("Error while reading file line by line:" 
	            + ex.getMessage());                      
	    }
	    
	    System.out.println("Graph read from file. " + num_nodes + " nodes, " + this.NumEdges() + " edges.");
	}

}
