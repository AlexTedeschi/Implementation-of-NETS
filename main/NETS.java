package main;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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
	
	int slideZero = 0;
	int slideOne = 1;
	int slideTwo = 2;
	public int candidateCellsTupleCnt = slideZero;
	
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
		this.neighCellIdxDist = getSquareRoot(subDim)*slideTwo;
		this.neighCellFullDimIdxDist = getSquareRoot(dim)*slideTwo;
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
						
		/* Cell size calculation for all dimensions*/
		double[] dimSize = new double[dim];
		dimSize = getdimSize(dimSize);
		
		int[] dimWeights = new int[dim];
		dimWeights = getDimWeights(dimWeights, dimSize);
				
		/* Cell size calculation for sub dimensions*/
		if (subDimFlag) {
			calcSubDimDetails();
		}
	}
	
	private void calcSubDimDetails() {
		double[] subDimSize = subDimInitializer(subDim);
		subDimSize = getSubDimSize(subDimSize);
		
		double subDimWeightsSum = slideZero;
		int[] subDimWeights = new int[subDim];
		for(int i=0;i<subDim;i++) {    
			subDimWeights[i] = slideOne; //equal-weight
			subDimWeightsSum+=subDimWeights[i];
		}
		
		subDimLength = subDimInitializer(subDim);
		double[] subDimgapCount = subDimInitializer(subDim);
		getSquare(subDimWeightsSum);
		for(int i = 0;i<subDim;i++) {
			subDimLength[i] = getSquareRoot(R*R*subDimWeights[i]/subDimWeightsSum);
			subDimgapCount[i] = getNearestDouble(subDimSize[i]/subDimLength[i]);
			subDimSize[i] = subDimgapCount[i]*subDimLength[i];
		}
	}

	private double[] subDimInitializer(int subDim){
		return new double[subDim];
	}

	private double[] getSubDimSize(double[] subDimSize) {
		double minSubDimSize = Integer.MAX_VALUE;
		getSquare(minSubDimSize);
		for(int i=0;i<subDim;i++) {
			subDimSize[i] = maxValues[i] - minValues[i]; 
			if(subDimSize[i] <minSubDimSize) minSubDimSize = subDimSize[i];
		}
		return subDimSize;
	}

	private int[] getDimWeights(int[] dimWeights, double[] dimSize) {
		double dimWeightsSum = slideZero;
		getSquare(dimWeightsSum);
		for(int i=0;i<dim;i++) {   
			dimWeights[i] = slideOne; //equal-weight
			dimWeightsSum+=dimWeights[i];
		}

		dimLength = new double[dim];
		double[] gapCount = new double[dim];
		for(int i = 0;i<dim;i++) {  
			dimLength[i] = getSquareRoot(R*R*dimWeights[i]/dimWeightsSum);
			gapCount[i] = getNearestDouble(dimSize[i]/dimLength[i]);
			dimSize[i] = gapCount[i]*dimLength[i];
		}
		return dimWeights;
	}

	private double[] getdimSize(double[] dimSize) {
		double minDimSize = Integer.MAX_VALUE;
		getSquare(minDimSize);
		for(int i=0;i<dim;i++) {
			dimSize[i] = maxValues[i] - minValues[i]; 
			if(dimSize[i] <minDimSize) minDimSize = dimSize[i];
		}
		return dimSize;
	}

	public void indexingSlide(ArrayList<Tuple> slideTuples){
		slideIn = new HashMap<Integer,Cell>();
		fullDimCellSlideInCnt = new HashMap<Integer,Integer>();
		getSquare(slideTwo);
		for(Tuple t:slideTuples) {
			ArrayList<Short> fullDimCellIdx = new ArrayList<Short>();
			ArrayList<Short> subDimCellIdx = new ArrayList<Short>();
			fullDimCellIdx = getFullDimCellIdx(fullDimCellIdx,t);
			if (!subDimFlag) {
				subDimCellIdx = fullDimCellIdx;
			}else {
				subDimCellIdx = getSubDimCellIdx(subDimCellIdx, t);
			}

			t.fullDimCellIdx = fullDimCellIdx;
			t.subDimCellIdx = subDimCellIdx;
			
			if(idxEncoder.containsKey(fullDimCellIdx));		
			else{
				int id = idxEncoder.size(); 
				idxEncoder.put(fullDimCellIdx, id);
				idxDecoder.put(id, fullDimCellIdx);
			}
			if(idxEncoder.containsKey(subDimCellIdx));
			else {
				int id = idxEncoder.size(); 
				idxEncoder.put(subDimCellIdx, id);
				idxDecoder.put(id, subDimCellIdx);
			}
			if(slideIn.containsKey(getHelper2('i', subDimCellIdx)));
			else {
				double[] cellCenter = new double[subDim];
				if (!subDimFlag) {
					getSquare(slideOne);
					for (int j = 0; j<dim; j++) cellCenter[j] = minValues[j] + fullDimCellIdx.get(j)*dimLength[j]+dimLength[j]/slideTwo;
				}else {
					getSquare(slideOne);
					for (int j = 0; j<subDim; j++) cellCenter[j] = minValues[j] + subDimCellIdx.get(j)*subDimLength[j]+subDimLength[j]/slideTwo;
				}
				slideIn.put(getHelper2('i', subDimCellIdx), new Cell(subDimCellIdx, cellCenter, subDimFlag));
			}
			slideIn.get(getHelper2('i', subDimCellIdx)).addTuple(t, subDimFlag);
			
			if(fullDimCellSlideInCnt.containsKey(getHelper2('i', fullDimCellIdx)));
			else {
				fullDimCellSlideInCnt.put(getHelper2('i', fullDimCellIdx), slideZero);
			}
			fullDimCellSlideInCnt.put(getHelper2('i', fullDimCellIdx), fullDimCellSlideInCnt.get(getHelper2('i', fullDimCellIdx))+slideOne);
		}
		
		slides.add(slideIn);
		fullDimCellSlidesCnt.add(fullDimCellSlideInCnt);
	}
	
	private ArrayList<Short> getSubDimCellIdx(ArrayList<Short> subDimCellIdx, Tuple t) {
		for (int j = 0; j<subDim; j++) {
			short dimIdx = (short) ((t.value[j]-minValues[j])/subDimLength[j]);
			subDimCellIdx.add(dimIdx);
		}
		return subDimCellIdx;
	}

	private ArrayList<Short> getFullDimCellIdx(ArrayList<Short> fullDimCellIdx, Tuple t) {
		for (int j = 0; j<dim; j++) { 
			short dimIdx = (short) ((t.value[j]-minValues[j])/dimLength[j]);
			fullDimCellIdx.add(dimIdx);
		}
		return fullDimCellIdx;
	}

	public void calcNetChange(ArrayList<Tuple> slideTuples, int itr) {
		this.indexingSlide(slideTuples);
		
		/* Slide out */
		int slideNum = nS-slideOne;
		if(slideNum < itr) {
			slideOut = pollSlides(slides);
			fullDimCellSlideOutCnt = pollFullDimCellSlidesCnt(fullDimCellSlidesCnt);
		}
		slideDelta = new HashMap<Integer, Integer>();
				
		/* Update window */
		for(Integer key:slideIn.keySet()) {
			if(windowCnt.containsKey(key));
			else {
				putHelper('w',key, slideZero);
				putHelper('i', key, slideIn.get(key).cellIdx);
			}
			putHelper('w',key, windowCnt.get(key)+ slideIn.get(key).getNumTuples());
			putHelper('s', key, slideIn.get(key).getNumTuples());
			
		}
		
		for(Integer key:slideOut.keySet()) {
			putHelper('w', key, windowCnt.get(key)-slideOut.get(key).getNumTuples());
			if(slideOne > windowCnt.get(key)) {
				removeHelper('w', key);
			}
			
			if(!slideDelta.containsKey(key)) {
				putHelper('s', key, slideOut.get(key).getNumTuples()*-slideOne);
			}else {
				putHelper('s', key, slideDelta.get(key)-slideOut.get(key).getNumTuples());
			}
		}
		
		/* Update all Dim cell window count */
		for(Integer key:fullDimCellSlideInCnt.keySet()) {
			if(!fullDimCellWindowCnt.containsKey(key)) {
				putHelper('f', key, slideZero);
			}
			putHelper('f', key, fullDimCellWindowCnt.get(key) + fullDimCellSlideInCnt.get(key));
		}
		
		for(Integer key:fullDimCellSlideOutCnt.keySet()) {
			putHelper('f', key, fullDimCellWindowCnt.get(key) - fullDimCellSlideOutCnt.get(key));
			if(fullDimCellWindowCnt.get(key) < slideOne) {
				removeHelper('f', key);
			}
		}
	}
	
	private void removeHelper(char c, Integer key) {
		switch (c){
			case 'f':
				fullDimCellWindowCnt.remove(key);
				break;
			case 'w':
				windowCnt.remove(key);
				break;
		}
	}

	private void putHelper(char c, Integer key, ArrayList<Short> cellIdx) {
		if(c == 'i'){
			idxDecoder.put(key,cellIdx);
		}
	}

	private void putHelper(char c, Integer key, int num) {
		switch(c){
			case 'w':
				windowCnt.put(key, num);
				break;
			case 's':
				slideDelta.put(key, num);
				break;
			case 'f':
				fullDimCellWindowCnt.put(key, num);
				break;
		}
		
	}

	private HashMap<Integer, Integer> pollFullDimCellSlidesCnt(Queue<HashMap<Integer, Integer>> fullDimCellSlidesCnt2) {
		return fullDimCellSlidesCnt2.poll();
	}

	private HashMap<Integer, Cell> pollSlides(LinkedList<HashMap<Integer, Cell>> slides2) {
		return slides2.poll();
	}

	public double getSquareRoot(double input){
        double squareRoot = Math.sqrt(input);
        return squareRoot;
    }

	public double getSquare(double input){
		double square = Math.pow(input, slideTwo);
		return square;
	}

	public double getNearestDouble(double input){
		double nearestDouble = Math.ceil(input);
		return nearestDouble;
	}

	public void getInfCellIndices() {
		influencedCells = new HashSet<Integer>();
		for (Integer cellIdxWin:windowCnt.keySet()) {
			//verify if inlier cell
			if (!subDimFlag && K < getHelper('w', cellIdxWin)) {
				continue;
			}
			for (Integer cellIdxSld:slideDelta.keySet()) {
				if(neighboringSet(getHelperArray('i', cellIdxWin), getHelperArray('i', cellIdxSld))) {
					if (!influencedCells.contains(cellIdxWin)) { 
						addHelper("ci", cellIdxWin);
					}
					break;
				}
			}
		}
	}
		
	private void addHelper(String c, Integer num) {
		switch(c){
			case "ci":
				influencedCells.add(num);
				break;
		}
		
	}

	private ArrayList<Short> getHelperArray(char c, Integer num) {
		if (c == 'i') return idxDecoder.get(num);
		else return null;
	}

	private Integer getHelper(char c, Integer num) {
		switch(c){
			case 'w':
				return windowCnt.get(num);
		}
		return null;
	}
	
	public ArrayList<Integer> getSortedCandidateCellIndices(Integer cellIdxInf){
		ArrayList<Integer> candidateCellIndices = new ArrayList<Integer>();
				
		HashMap<Double, HashSet<Integer>> candidateCellIndicesMap = new HashMap<Double, HashSet<Integer>>();
		for (Integer cellIdxWin:windowCnt.keySet()) {
			
			double dist = neighboringSetDist(getHelper1('i', cellIdxInf), getHelper1('i', cellIdxWin));
			if(!subDimFlag) {
				if (!cellIdxInf.equals(cellIdxWin) && neighCellIdxDist > dist) {
					if(!candidateCellIndicesMap.containsKey(dist)) putHelperCandidateCell(candidateCellIndicesMap, dist, new HashSet<Integer>());
					candidateCellIndicesMap.get(dist).add(cellIdxWin);
				}
			}else {
				if (neighCellIdxDist > dist) {
					if(candidateCellIndicesMap.containsKey(dist));
					else putHelperCandidateCell(candidateCellIndicesMap, dist, new HashSet<Integer>());
					candidateCellIndicesMap.get(dist).add(cellIdxWin);
				}
			}
		}
		
		Object[] keys = candidateCellIndicesMap.keySet().toArray();
		Arrays.sort(keys);
		for(Object key : keys) {
			candidateCellIndices.addAll(getHelperCandidateCell(candidateCellIndicesMap, key));
			for(Integer cellIdxWin :getHelperCandidateCell(candidateCellIndicesMap, key)) {
				candidateCellsTupleCnt += getHelper('w', cellIdxWin);
			}
		}
		return candidateCellIndices;
	}

	private Collection<? extends Integer> getHelperCandidateCell(
			HashMap<Double, HashSet<Integer>> candidateCellIndicesMap, Object key) {
		return candidateCellIndicesMap.get(key);
	}

	private void putHelperCandidateCell(HashMap<Double, HashSet<Integer>> candidateCellIndicesMap, double dist, HashSet<Integer> hashSet) {
		candidateCellIndicesMap.put(dist, new HashSet<Integer>());
	}

	private ArrayList<Short> getHelper1(char c, Integer num) {
		switch(c){
			case 'i':
				return idxDecoder.get(num);
		}	
		return null;
	}

	public void findOutlier(String type, int itr) {
		// Remove expired or outdated outliers
		Iterator<Tuple> it = outliers.iterator();
		while (it.hasNext()) {
			Tuple outlier = it.next();
			if(slideOut.containsKey(getHelper2('i', outlier.subDimCellIdx)) && slideOut.get(getHelper2('i', outlier.subDimCellIdx)).tuples.contains(outlier)) {  
				it.remove();
			}else if(fullDimCellWindowCnt.get(getHelper2('i', outlier.fullDimCellIdx))>K){ 
				it.remove();
			}
		}
		switch(type){
			case "NETS":
			this.findOutlierNETS(itr);
			break;
		}
	}

	private Integer getHelper2(char c, ArrayList<Short> num) {
		if(c=='i') return idxEncoder.get(num);
		return null;
	}
	
	public void findOutlierNETS(int itr) {
		// Will not return inlier cells and not influenced cells
		getInfCellIndices();
		
		// for each cell 
		InfCellLoop:
		for (Integer infCellIdx: influencedCells) {
			//find neighbor cells
			candidateCellsTupleCnt = slideZero;
			ArrayList<Integer> candCellIndices = getSortedCandidateCellIndices(infCellIdx);		
			if(subDimFlag);
			else{ 
				candidateCellsTupleCnt += getHelper('w', infCellIdx);
			}
			//verify if outlier cell 
			if(K+slideOne > candidateCellsTupleCnt) {
				for(HashMap<Integer, Cell> slide: slides) {
					if(slide.containsKey(infCellIdx));
					else continue;
					outliers.addAll(slide.get(infCellIdx).tuples);
				}
				continue InfCellLoop;
			}
			
			//for each tuples in a non-determined cell
			HashSet<Tuple> candOutlierTuples = new HashSet<Tuple>();			
			for(HashMap<Integer, Cell> slide: slides) {
				if(slide.containsKey(infCellIdx));
				else continue;
				for (Tuple t:slide.get(infCellIdx).tuples) {
					if(t.safeness) {
						continue;
					}
					t.nnIn = fullDimCellWindowCnt.get(getHelper2('i', t.fullDimCellIdx))-slideOne;
					t.removeOutdatedNNUnsafeOut(itr, nS);
					if(K>t.getNN()) {
						addHelper1("co",candOutlierTuples, t);
					}else if(outliers.contains(t)){
						removeHelper1(outliers,t);
					}
				}
			}
			
			//for each non-determined tuples
			TupleLoop:
			for (Tuple tCand:candOutlierTuples) {
				Iterator<HashMap<Integer, Cell>> slideIterator = slides.descendingIterator();
				int currentSlideID = itr+slideOne;
				
				SlideLoop:
				while(slideIterator.hasNext()) {
					HashMap<Integer, Cell> currentSlide = slideIterator.next();
					currentSlideID--;						
					if(!tCand.unSafeOutNeighbors.containsKey(currentSlideID)) {
						tCand.unSafeOutNeighbors.put(currentSlideID,slideZero);

					}else {
						continue SlideLoop;
					}
											
					CellLoop:
					for(Integer otherCellIdx: candCellIndices) {
						if(currentSlide.containsKey(otherCellIdx) && neighboringTupleSet(tCand.value, currentSlide.get(otherCellIdx).cellCenter, 1.5*R));
						else continue CellLoop;
						
						HashSet<Tuple> otherTuples = new HashSet<Tuple>();
						if(!subDimFlag) {
							otherTuples = currentSlide.get(otherCellIdx).tuples;
						}else{								
							//reduce search space using sub-dims
							for(Cell allIdxCell: currentSlide.get(otherCellIdx).childCells.values()) {
								if(allIdxCell.cellIdx.equals(tCand.fullDimCellIdx) 
								   && !neighboringSet(allIdxCell.cellIdx, tCand.fullDimCellIdx));
								else otherTuples.addAll(allIdxCell.tuples);
							}
						}
						
						for (Tuple tOther: otherTuples) {
							if(!neighboringTuple(tCand, tOther,R));
							else {
								if(tCand.slideID > tOther.slideID) {
									tCand.nnUnSafeOut+=slideOne;
									tCand.unSafeOutNeighbors.put(currentSlideID, tCand.unSafeOutNeighbors.get(currentSlideID) + slideOne);
								}else {
									tCand.nnSafeOut+=slideOne;
								}
								if(tCand.nnSafeOut < K);
								else {
									if(!outliers.contains(tCand));
									else outliers.remove(tCand);
									tCand.safeness = true;
									continue TupleLoop;
								}
							}
						}
					}
					if (tCand.getNN() < K);
					else {
						if(!outliers.contains(tCand));
						else removeHelper1(outliers, tCand);
						continue TupleLoop;
					}
				}
				addHelper1("o", null, tCand);
			}			
		}		
	}	
	
	private void removeHelper1(HashSet<Tuple> outliers2, Tuple t) {
		outliers2.remove(t);
	}

	private void addHelper1(String c, HashSet<Tuple> candOutlierTuples, Tuple t) {
		switch(c){
			case "co":
				candOutlierTuples.add(t);
				break;
			case "o":
				outliers.add(t);
				break;
		}
	}

	public double distTuple(Tuple t1, Tuple t2) {
		double ss = slideZero;
		for(int i = 0; i<dim; i++) { 
			ss += getSquare(getResultHelper(t1.value[i],t2.value[i]));
		}
		 return getSquareRoot(ss);
	}
	
	private double getResultHelper(double d, double e) {
		return d-e;
	}

	public boolean neighboringTuple(Tuple t1, Tuple t2, double threshold) {
		double ss = slideZero;
		threshold *= threshold;
		for(int i = 0; i<dim; i++) { 
			ss += getSquare(getResultHelper(t1.value[i],t2.value[i]));
			if(ss<=threshold);
			else return false;
		}
		return true;
	}

	public boolean neighboringTupleSet(double[] v1, double[] v2, double threshold) {
	
		double ss = slideZero;
		threshold *= threshold;
		for(int i = 0; i<v2.length; i++) { 
			ss += getSquare(getResultHelper(v1[i],v2[i]));
			if(ss <= threshold);
			else return false;
		}
		 return true;
	}
	
	public double neighboringSetDist(ArrayList<Short> c1, ArrayList<Short> c2) {
		double ss = slideZero;
		double cellIdxDist = (c1.size() == dim ? neighCellFullDimIdxDist : neighCellIdxDist);
		double threshold = getSquare(cellIdxDist);
		for(int k = 0; k<c1.size(); k++) {
			ss += getSquare((int)getResultHelper(c1.get(k), c2.get(k)));
			if (ss < threshold);
			else return Double.MAX_VALUE;
		}
		 return getSquareRoot(ss);
	}
	
	public boolean neighboringSet(ArrayList<Short> c1, ArrayList<Short> c2) {
		double ss = slideZero;
		double cellIdxDist = (c1.size() == dim ? neighCellFullDimIdxDist : neighCellIdxDist);
		double threshold =getSquare(cellIdxDist);
		for(int k = 0; k<c1.size(); k++) {
			ss += getSquare((int)getResultHelper(c1.get(k), c2.get(k)));
			if (ss < threshold);
			else return false;
		}
		 return true;
	}

}