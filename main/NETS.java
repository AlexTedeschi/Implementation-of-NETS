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
	public double radius;
	public int k_neighbors;
	public int dimension;
	public int sub_dimension;
	public boolean subdimension_flag;
	public int slide_size;
	public int window_size;
	public int number_of_slides;
	public int number_of_windows;
	public double neighCellIdxDist;
	public double neighCellFullDimIdxDist;
	public double[] maxValues;
	public double[] minValues;
	public double[] dimension_Length;
	public double[] subDimension_Length;
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
	
	public NETS(int dimension, int sub_dimension, double radius, int K, int slide_size, int window_size, int number_of_windows, double[] maxValues, double[] minValues) {
		this.dimension = dimension;
		this.sub_dimension = sub_dimension;
		this.subdimension_flag = dimension != sub_dimension;
		this.radius = radius;
		this.k_neighbors = K;
		this.slide_size = slide_size;
		this.window_size = window_size;
		this.number_of_windows = number_of_windows;
		this.number_of_slides = window_size / slide_size;
		this.neighCellIdxDist = getSquareRoot(sub_dimension)*slideTwo;
		this.neighCellFullDimIdxDist = getSquareRoot(dimension)*slideTwo;
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
		double[] dimSize = new double[dimension];
		dimSize = getdimSize(dimSize);
		
		int[] dimWeights = new int[dimension];
		dimWeights = getDimWeights(dimWeights, dimSize);
				
		/* Cell size calculation for sub dimensions*/
		if (subdimension_flag) {
			calcSubDimDetails();
		}
	}
	
	private void calcSubDimDetails() {
		double[] subDimSize = subDimInitializer(sub_dimension);
		subDimSize = getSubDimSize(subDimSize);
		
		double subDimWeightsSum = slideZero;
		int[] subDimWeights = new int[sub_dimension];
		for(int i = 0; i< sub_dimension; i++) {
			subDimWeights[i] = slideOne; //equal-weight
			subDimWeightsSum+=subDimWeights[i];
		}
		
		subDimension_Length = subDimInitializer(sub_dimension);
		double[] subDimgapCount = subDimInitializer(sub_dimension);
		getSquare(subDimWeightsSum);
		for(int i = 0; i< sub_dimension; i++) {
			subDimension_Length[i] = getSquareRoot(radius * radius *subDimWeights[i]/subDimWeightsSum);
			subDimgapCount[i] = getNearestDouble(subDimSize[i]/subDimension_Length[i]);
			subDimSize[i] = subDimgapCount[i]*subDimension_Length[i];
		}
	}

	private double[] subDimInitializer(int subDim){
		return new double[subDim];
	}

	private double[] getSubDimSize(double[] subDimSize) {
		double minSubDimSize = Integer.MAX_VALUE;
		getSquare(minSubDimSize);
		for(int i = 0; i< sub_dimension; i++) {
			subDimSize[i] = maxValues[i] - minValues[i]; 
			if(subDimSize[i] <minSubDimSize) minSubDimSize = subDimSize[i];
		}
		return subDimSize;
	}

	private int[] getDimWeights(int[] dimWeights, double[] dimSize) {
		double dimWeightsSum = slideZero;
		getSquare(dimWeightsSum);
		for(int i = 0; i< dimension; i++) {
			dimWeights[i] = slideOne; //equal-weight
			dimWeightsSum+=dimWeights[i];
		}

		dimension_Length = new double[dimension];
		double[] gapCount = new double[dimension];
		for(int i = 0; i< dimension; i++) {
			dimension_Length[i] = getSquareRoot(radius * radius *dimWeights[i]/dimWeightsSum);
			gapCount[i] = getNearestDouble(dimSize[i]/ dimension_Length[i]);
			dimSize[i] = gapCount[i]* dimension_Length[i];
		}
		return dimWeights;
	}

	private double[] getdimSize(double[] dimSize) {
		double minDimSize = Integer.MAX_VALUE;
		getSquare(minDimSize);
		for(int i = 0; i< dimension; i++) {
			dimSize[i] = maxValues[i] - minValues[i]; 
			if(dimSize[i] <minDimSize) minDimSize = dimSize[i];
		}
		return dimSize;
	}

	// Cell level outlier detection
	public void indexingSlide(ArrayList<Tuple> slideTuples){
		slideIn = new HashMap<Integer,Cell>();
		fullDimCellSlideInCnt = new HashMap<Integer,Integer>();
		getSquare(slideTwo);
		for(Tuple t:slideTuples) {
			ArrayList<Short> fullDimCellIdx = new ArrayList<Short>();
			ArrayList<Short> subDimCellIdx = new ArrayList<Short>();
			fullDimCellIdx = getFullDimCellIdx(fullDimCellIdx,t);
			if (!subdimension_flag) {
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
				double[] cellCenter = new double[sub_dimension];
				if (!subdimension_flag) {
					getSquare(slideOne);
					for (int j = 0; j< dimension; j++) cellCenter[j] = minValues[j] + fullDimCellIdx.get(j)* dimension_Length[j]+ dimension_Length[j]/slideTwo;
				}else {
					getSquare(slideOne);
					for (int j = 0; j< sub_dimension; j++) cellCenter[j] = minValues[j] + subDimCellIdx.get(j)*subDimension_Length[j]+subDimension_Length[j]/slideTwo;
				}
				slideIn.put(getHelper2('i', subDimCellIdx), new Cell(subDimCellIdx, cellCenter, subdimension_flag));
			}
			slideIn.get(getHelper2('i', subDimCellIdx)).addTuple(t, subdimension_flag);
			
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
		for (int j = 0; j< sub_dimension; j++) {
			short dimIdx = (short) ((t.value[j]-minValues[j])/subDimension_Length[j]);
			subDimCellIdx.add(dimIdx);
		}
		return subDimCellIdx;
	}

	private ArrayList<Short> getFullDimCellIdx(ArrayList<Short> fullDimCellIdx, Tuple t) {
		for (int j = 0; j< dimension; j++) {
			short dimIdx = (short) ((t.value[j]-minValues[j])/ dimension_Length[j]);
			fullDimCellIdx.add(dimIdx);
		}
		return fullDimCellIdx;
	}

	// For calculating slide net effect change (Preprocessing)
	public void calcNetChange(ArrayList<Tuple> slideTuples, int itr) {
		this.indexingSlide(slideTuples);
		
		/* Slide out */
		int slideNum = number_of_slides -slideOne;
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
			if (!subdimension_flag && k_neighbors < getHelper('w', cellIdxWin)) {
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
			if(!subdimension_flag) {
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

	// Find outliers in the data stream
	public void findOutlier(String type, int itr) {
		// Remove expired or outdated outliers
		Iterator<Tuple> it = outliers.iterator();
		while (it.hasNext()) {
			Tuple outlier = it.next();
			if(slideOut.containsKey(getHelper2('i', outlier.subDimCellIdx)) && slideOut.get(getHelper2('i', outlier.subDimCellIdx)).tuples.contains(outlier)) {  
				it.remove();
			}else if(fullDimCellWindowCnt.get(getHelper2('i', outlier.fullDimCellIdx))>k_neighbors){
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

	// This function is used to find outliers in a given data stream
	public void findOutlierNETS(int itr) {
		// Will not return inlier cells and not influenced cells
		getInfCellIndices();

		//Point level outlier detection
		InfCellLoop:
		for (Integer infCellIdx: influencedCells) {
			//find neighbor cells
			candidateCellsTupleCnt = slideZero;
			ArrayList<Integer> candCellIndices = getSortedCandidateCellIndices(infCellIdx);		
			if(subdimension_flag);
			else{ 
				candidateCellsTupleCnt += getHelper('w', infCellIdx);
			}
			//verify if outlier cell 
			if(k_neighbors+slideOne > candidateCellsTupleCnt) {
				for(HashMap<Integer, Cell> slide: slides) {
					if(slide.containsKey(infCellIdx));
					else continue;
					outliers.addAll(slide.get(infCellIdx).tuples);
				}
				continue InfCellLoop;
			}
			

			HashSet<Tuple> candOutlierTuples = new HashSet<Tuple>();			
			for(HashMap<Integer, Cell> slide: slides) {
				if(slide.containsKey(infCellIdx));
				else continue;
				for (Tuple t:slide.get(infCellIdx).tuples) {
					if(t.safeness) {
						continue;
					}
					t.nnIn = fullDimCellWindowCnt.get(getHelper2('i', t.fullDimCellIdx))-slideOne;
					t.removeOutdatedNNUnsafeOut(itr, number_of_slides);
					if(k_neighbors>t.getNN()) {
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
						if(currentSlide.containsKey(otherCellIdx) && neighboringTupleSet(tCand.value, currentSlide.get(otherCellIdx).cellCenter, 1.5* radius));
						else continue CellLoop;
						
						HashSet<Tuple> otherTuples = new HashSet<Tuple>();
						if(!subdimension_flag) {
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
							if(!neighboringTuple(tCand, tOther, radius));
							else {
								if(tCand.slideID > tOther.slideID) {
									tCand.nnUnSafeOut+=slideOne;
									tCand.unSafeOutNeighbors.put(currentSlideID, tCand.unSafeOutNeighbors.get(currentSlideID) + slideOne);
								}else {
									tCand.nnSafeOut+=slideOne;
								}
								if(tCand.nnSafeOut < k_neighbors);
								else {
									if(!outliers.contains(tCand));
									else outliers.remove(tCand);
									tCand.safeness = true;
									continue TupleLoop;
								}
							}
						}
					}
					if (tCand.getNN() < k_neighbors);
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
		for(int i = 0; i< dimension; i++) {
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
		for(int i = 0; i< dimension; i++) {
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
		double cellIdxDist = (c1.size() == dimension ? neighCellFullDimIdxDist : neighCellIdxDist);
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
		double cellIdxDist = (c1.size() == dimension ? neighCellFullDimIdxDist : neighCellIdxDist);
		double threshold =getSquare(cellIdxDist);
		for(int k = 0; k<c1.size(); k++) {
			ss += getSquare((int)getResultHelper(c1.get(k), c2.get(k)));
			if (ss < threshold);
			else return false;
		}
		 return true;
	}

}