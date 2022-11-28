package src.main;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;

public class NETS {
    public double R;
	public int K;
	public int dim;
	public int subDim;
	public boolean subDimFlag;
	public int S;
	public int W;
	public int nS;
	public int nW;
	public double neighCellIdxDist;
	public double neighCellFullDimIdxDist;
	public double[] maxValues;
	public double[] minValues;
	public double[] dimLength;
	public double[] subDimLength;
	
	public HashMap<Integer,Cell> slideIn;
	public HashMap<Integer,Cell> slideOut;
	public HashMap<Integer,Integer> windowCnt;
	public HashMap<Integer,ArrayList<Short>> idxDecoder;
	public HashMap<ArrayList<Short>,Integer> idxEncoder;
	public HashMap<Integer,Integer> slideDelta;
	public HashSet<Tuple> outliers;
	public HashMap<Integer,Integer> fullDimCellWindowCnt;
	public LinkedList<HashMap<Integer,Cell>> slides; 
	public HashMap<Integer,Integer> fullDimCellSlideInCnt;
	public HashMap<Integer,Integer> fullDimCellSlideOutCnt;
	public Queue<HashMap<Integer,Integer>> fullDimCellSlidesCnt; 
	public HashSet<Integer> influencedCells;
	
	public int candidateCellsTupleCnt = 0;

    public NETS(int dim, int subDim, double R, int K, int S, int W, int nW, double[] maxValues, double[] minValues) {
		this.dim = dim;
		this.subDim = subDim;
		this.subDimFlag = dim != subDim;
		this.R = R;
		this.K = K;
		this.S = S;
		this.W = W;
		this.nW = nW;
		this.nS = W/S;
		// this.neighCellIdxDist = Math.sqrt(subDim)*2;
		this.neighCellIdxDist = getSquareRoot(subDim)*2;
		// this.neighCellFullDimIdxDist = Math.sqrt(dim)*2;
		this.neighCellFullDimIdxDist = getSquareRoot(dim)*2;
		this.maxValues = maxValues;
		this.minValues = minValues;
		
		this.windowCnt = new HashMap<Integer,Integer>();
		this.slides = new LinkedList<HashMap<Integer,Cell>>();
		this.slideOut = new HashMap<Integer,Cell>();
		this.idxDecoder = new HashMap<Integer, ArrayList<Short>>();
		this.idxEncoder = new HashMap<ArrayList<Short>,Integer>();		
		this.fullDimCellWindowCnt = new HashMap<Integer,Integer>();
		this.fullDimCellSlidesCnt = new LinkedList<HashMap<Integer,Integer>>();
		this.fullDimCellSlideOutCnt = new HashMap<Integer,Integer>();
				
		this.outliers = new HashSet<Tuple>();
						
		/* Cell size calculation for all dim*/
		double minDimSize = Integer.MAX_VALUE;
		double[] dimSize = new double[dim];
		for(int i=0;i<dim;i++) {
			dimSize[i] = maxValues[i] - minValues[i]; 
			if(dimSize[i] <minDimSize) minDimSize = dimSize[i];
		}
		
		double dimWeightsSum = 0;
		int[] dimWeights = new int[dim];
		for(int i=0;i<dim;i++) {   
			//dimWeights[i] = dimSize[i]/minDimSize; //relative-weight
			dimWeights[i] = 1; //equal-weight
			dimWeightsSum+=dimWeights[i];
		}
		
		dimLength = new double[dim];
		double[] gapCount = new double[dim];
		for(int i = 0;i<dim;i++) {  
			// dimLength[i] = Math.sqrt(R*R*dimWeights[i]/dimWeightsSum);
			dimLength[i] = getSquareRoot(R*R*dimWeights[i]/dimWeightsSum);
			gapCount[i] = Math.ceil(dimSize[i]/dimLength[i]);
			dimSize[i] = gapCount[i]*dimLength[i];
		}
		
		/* Cell size calculation for sub dim*/
		if (subDimFlag) {
			double minSubDimSize = Integer.MAX_VALUE;
			double[] subDimSize = new double[subDim];
			for(int i=0;i<subDim;i++) {
				subDimSize[i] = maxValues[i] - minValues[i]; 
				if(subDimSize[i] <minSubDimSize) minSubDimSize = subDimSize[i];
			}
			
			double subDimWeightsSum = 0;
			int[] subDimWeights = new int[subDim];
			for(int i=0;i<subDim;i++) {    
				//subDimWeights[i] = subDimSize[i]/minSubDimSize; //relative-weight
				subDimWeights[i] = 1; //equal-weight
				subDimWeightsSum+=subDimWeights[i];
			}
			
			subDimLength = new double[subDim];
			double[] subDimgapCount = new double[subDim];
			for(int i = 0;i<subDim;i++) {   
				// subDimLength[i] = Math.sqrt(R*R*subDimWeights[i]/subDimWeightsSum);
				subDimLength[i] = getSquareRoot(R*R*subDimWeights[i]/subDimWeightsSum);
				subDimgapCount[i] = Math.ceil(subDimSize[i]/subDimLength[i]);
				subDimSize[i] = subDimgapCount[i]*subDimLength[i];
			}
		}

	}

    public double getSquareRoot(double input){
        double squareRoot = Math.sqrt(input);
        return squareRoot;
    }
    
}
