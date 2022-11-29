package utils;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import main.Tuple;

public class StreamGenerator {
	private double[] maxValues;
	private double[] minValues;
	private BufferedReader br; 
	private String filePath;
	private List<Integer> priorityList;
	
	public StreamGenerator (String dataset, int random) throws FileNotFoundException {
		/* Datasets */
		switch (dataset) {
			default: case "TAO":
				filePath = "datasets/TAO.csv";
				this.maxValues = new double[]{75.39,101.68,30.191};
				this.minValues = new double[]{-9.99,-9.99,-9.999};
				this.priorityList = Arrays.asList(new Integer[] {1,2,0});
				sortPriority(random);
				break;
			case "STK":
				filePath = "datasets/Stock.csv";
				this.maxValues = new double[]{9930};
				this.minValues = new double[]{0};	
				this.priorityList = Arrays.asList(new Integer[] {0});
				sortPriority(random);
				break;
			case "GAU":
				filePath = "datasets/Gauss.csv";
				this.maxValues = new double[]{100.81};
				this.minValues = new double[]{-3.5042};
				this.priorityList = Arrays.asList(new Integer[] {0});
				sortPriority(random);
				break;
		}
	}
	
	public void sortPriority(int random){
		if(random>0) Collections.shuffle(this.priorityList);
		double[] new_maxValues = new double[priorityList.size()];
		double[] new_minValues = new double[priorityList.size()];
		for(int i=0; i< priorityList.size(); i++) {
			new_maxValues[i] = maxValues[this.priorityList.get(i)];
			new_minValues[i] = minValues[this.priorityList.get(i)];
		}
		this.maxValues = new_maxValues;
		this.minValues = new_minValues;
	}
	
	public ArrayList<Tuple> getNewSlideTuples(int itr, int S) throws IOException {
		ArrayList<Tuple> newSlide = new ArrayList<Tuple>();
		this.br = new BufferedReader(new FileReader(filePath));
		String line = br.readLine();
		int tid = 0;
		
		while(line!=null) {
			if(tid>=itr*S) {
				String[] rawValues = line.split(",");
				double[] value = new double[rawValues.length];
				
				int j = 0;
				for(int i: priorityList) {
					value[j] = Double.parseDouble(rawValues[i]);
					j++;
				}
				Tuple tuple = new Tuple(tid, itr, value);
				newSlide.add(tuple);
			}
			tid++;
			if(tid==(itr+1)*S) break;
			line = br.readLine();
		}
		return newSlide;
	}

	public double[] getMaxValues() {
		return this.maxValues;
	}
	public double[] getMinValues(){
		return this.minValues;
	}
	
	public void setPriorityList(Integer[] list){
		this.priorityList = Arrays.asList(list);
	}
	
}

