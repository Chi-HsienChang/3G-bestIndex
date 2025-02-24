/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include <iostream>
#include "chromosome.h"
#include "dsmga2.h"
#include "fastcounting.h"
#include "statistics.h"

#include <iomanip>
using namespace std;


DSMGA2::DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff) {


    previousFitnessMean = -INF;
    ell = n_ell;
    nCurrent = (n_nInitial/2)*2;  // has to be even
    
    isIncrease_gen01 = vector<double>(ell, 0);
    isIncrease_gen12 = vector<double>(ell, 0);


    ratio_1 = vector<vector<double>>(3, vector<double>(ell, 0));
    ratio_1_check = vector<vector<double>>(3, vector<double>(ell, 0));

    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0;
    Chromosome::lsnfe = 0;
    Chromosome::hitnfe = 0;
    Chromosome::hit = false;

    selectionPressure = 2;
    maxGen = n_maxGen;
    maxFe = n_maxFe;

    graph.init(ell);
    graph_size.init(ell);
     
    bestIndex = -1;
    masks = new list<int>[ell];
    selectionIndex = new int[nCurrent];
    orderN = new int[nCurrent];
    orderELL = new int[ell];
    population = new Chromosome[nCurrent];
    temp_population = new Chromosome[nCurrent+2];
    copy_population = new Chromosome[nCurrent+2];
    fastCounting = new FastCounting[ell];

    for (int i = 0; i < ell; i++)
        fastCounting[i].init(nCurrent);


    pHash.clear();
    for (int i=0; i<nCurrent; ++i) {
        population[i].initR(ell);
        double f = population[i].getFitness();
        pHash[population[i].getKey()] = f;
    }

    if (GHC) {
        for (int i=0; i < nCurrent; i++)
            population[i].GHC();
    }
}


DSMGA2::~DSMGA2 () {
    delete []masks;
    delete []orderN;
    delete []orderELL;
    delete []selectionIndex;
    delete []population;
    delete []fastCounting;
}



bool DSMGA2::isSteadyState () {

    if (stFitness.getNumber () <= 0)
        return false;

    if (previousFitnessMean < stFitness.getMean ()) {
        previousFitnessMean = stFitness.getMean () + EPSILON;
        return false;
    }

    return true;
}

bool DSMGA2::converged() {
    if (stFitness.getMax() == lastMax &&
        stFitness.getMean() == lastMean &&
        stFitness.getMin() == lastMin)
        convergeCount++;
    else
        convergeCount = 0;

    lastMax = stFitness.getMax();
    lastMean = stFitness.getMean();
    lastMin = stFitness.getMin();
    return (convergeCount > 300) ? true : false;
}

int DSMGA2::doIt (bool output) {
    generation = 0;
    while (!shouldTerminate ()) {
        oneRun (output);
       #ifdef DEBUG
       cin.get();
        #endif
    }
    return generation;
}

    //  ONEMAX:  0
    //  MK    :  1
    //  FTRAP :  2
    //  CYC   :  3
    //  NK    :  4
    //  SPIN  :  5
    //  SAT   :  6
    //  mixTRAP:  7
    //  THREEOPTIMUM:  8



void DSMGA2::oneRun (bool output) {


    if (generation <= 2) {
        if (CACHE)
            Chromosome::cache.clear();

        for (int i = 0; i < ell; ++i) {
            double count1 = 0;
            for (int j = 0; j < nCurrent; ++j) {
                if (population[j].getVal(i) == 1)
                    count1++;
            }
            ratio_1[generation][i] = count1 / nCurrent;
        }


        if (generation == 1)
        {
            for (int i = 0; i < ell; ++i) {  
                isIncrease_gen01[i] = ratio_1[1][i] - ratio_1[0][i];
            }
        }


        if (generation == 2) {
            for (int i = 0; i < ell; ++i) {
                isIncrease_gen12[i] = ratio_1[2][i] - ratio_1[1][i];     
            }

            double change01_12 = 0;

            for (int i = 0; i < ell; ++i) {
                if (isIncrease_gen01[i] > 0 && isIncrease_gen12[i] < 0 || isIncrease_gen01[i] < 0 && isIncrease_gen12[i] > 0) {
                    change01_12++;
                }        
            }

            if (change01_12 >= ell*0.01) {
                cluster = true;
            } else {
                cluster = false;
            }
                
        }

        mixing();   

    }else{
        
        if (!cluster)
        {
            mixing();  
        }{
            
            cout << "Generation: " << generation << endl;
            for (int i = 0; i < nCurrent; ++i) {
                copy_population[i] = population[i];
            }

            // cout << "population before cluster:" << endl;
            // for (int i = 0; i < nCurrent; ++i) {
            //     cout << i << "= "; 
            //     for (int j = 0; j < ell; ++j) {
            //        cout << population[i].getVal(j);
            //     }
            //     cout << endl;
            // }

            // cout << "copy_population before cluster:" << endl;
            // for (int i = 0; i < nCurrent; ++i) {
            //     cout << i << "= "; 
            //     for (int j = 0; j < ell; ++j) {
            //        cout << copy_population[i].getVal(j);
            //     }
            //     cout << endl;
            // }
        
            int original_nCurrent = nCurrent;

            std::vector<int> group1, group2, group3;
            std::vector<int> seeds(3, -1);

            seeds[0] = bestIndex;

            int maxdistance = -1;
            int hamm = -1;
            for(int i = 0; i < nCurrent; i++){
                if (i == bestIndex) continue;

                hamm = 0;
                for (int j = 0; j < ell; ++j) {
                    if (population[bestIndex].getVal(j) != population[i].getVal(j))
                        ++hamm;
                }

                if (hamm > maxdistance){
                    seeds[1] = i;
                    maxdistance = hamm;
                }
            }


            int mindistance = ell;
            int hamm2 = -1;
            for(int i = 0; i < nCurrent; i++){
                if (i == bestIndex) continue;

                hamm2 = 0;
                for (int j = 0; j < ell; ++j) {
                    if (population[bestIndex].getVal(j) != population[i].getVal(j))
                        ++hamm2;
                }

                if (abs(hamm2 - (maxdistance/2)) < mindistance){
                    seeds[2] = i;
                    mindistance = abs(hamm2 - (maxdistance/2));
                    // cout << "mindistance: " << mindistance << endl;
                }
            }  

            vector<vector<int>> groups(3);
            groups[0].push_back(seeds[0]);
            groups[1].push_back(seeds[1]);
            groups[2].push_back(seeds[2]);

            for (size_t i = 0; i < nCurrent; ++i) {
                if (i == seeds[0] || i == seeds[1] || i == seeds[2]) continue; // 跳過種子

                int bestGroup = 0;
                int bestWeightedSimilarity = -1;
                for (int g = 0; g < 3; ++g) {
                    double groupSimilarity = 0;

                    // 計算該群組的總相似性
                    for (int idx : groups[g]) {
                        int sim = 0;
                        for (size_t geneIdx = 0; geneIdx < ell; ++geneIdx) {
                            if (population[i].getVal(geneIdx) == population[idx].getVal(geneIdx)) {
                                sim++;
                            }
                        }
                        groupSimilarity += sim;
                    }

                
                    double weightedSimilarity = groups[g].empty() ? groupSimilarity : (groupSimilarity / groups[g].size());

                
                    if (weightedSimilarity > bestWeightedSimilarity) {
                        bestWeightedSimilarity = weightedSimilarity;
                        bestGroup = g;
                    }
                }
                
                groups[bestGroup].push_back(i);
            }
            
            // cout << "groups[0].size: "<< groups[0].size() <<", groups[1].size: "<< groups[1].size() <<", groups[2].size: "<< groups[2].size() << endl;
            // for (int idx : groups[0]) {
            //     cout << idx << " ";
            // }
            // cout << endl;

            // for (int i = 0; i < groups[0].size(); ++i) {
            //     for (int j = 0; j < ell; j++) {
            //         cout << population[groups[0][i]].getVal(j);
            //     }
            //     cout << endl;        
            // }
            // cout << "---------------" << endl;

            // cout << "groups[1].size: "<< groups[1].size() << endl;
            // for (int idx : groups[1]) {
            //     cout << idx << " ";
            // }
            // cout << endl;

            // for (int i = 0; i < groups[1].size(); ++i) {
            //     for (int j = 0; j < ell; j++) {
            //         cout << population[groups[1][i]].getVal(j);
            //     }
            //     cout << endl;        
            // }

            // cout << "---------------" << endl;

            // cout << "groups[2].size: "<< groups[2].size() << endl;
            // for (int idx : groups[2]) {
            //     cout << idx << " ";
            // }
            // cout << endl;

            // for (int i = 0; i < groups[2].size(); ++i) {
            //     for (int j = 0; j < ell; j++) {
            //         cout << population[groups[2][i]].getVal(j);
            //     }
            //     cout << endl;        
            // }

            bool copy_population_active = true;

            int i_count = 0;
            for (int G = 0; G < 3; ++G) {
                nCurrent = groups[G].size();


                if (copy_population_active){
                    for (int i = 0; i < nCurrent; ++i) {
                        population[i] = copy_population[groups[G][i]]; //修改後的寫法
                    }
                }else{
                     for (int i = 0; i < nCurrent; ++i) {
                        population[i] = population[groups[G][i]]; //原本的寫法
                    }                   
                }
                


                // cout << "G = " << G << " nCurrent = " << nCurrent << endl;
                // for (int i = 0; i < nCurrent; ++i) {
                //     for (int j = 0; j < ell; j++) {
                //         cout << population[i].getVal(j);
                //     }
                //     cout << endl;
                // }                

                if (CACHE)
                    Chromosome::cache.clear();

                mixing();

            
                for (int i = 0; i < nCurrent; ++i) {
                    if (i+i_count >= original_nCurrent) break;
                    temp_population[i+i_count] = population[i];
                }

                i_count += nCurrent;

                if (G == 0) {
                    nCurrentA = nCurrent;
                } else if (G == 1) {
                    nCurrentB = nCurrent;
                } else {
                    nCurrentC = nCurrent;
                }

            }

            nCurrent = original_nCurrent;
            // 將群組中的染色體放回 population
            for (int i = 0; i < nCurrent; ++i) {
                population[i] = temp_population[i];
            }
        }
    }


    // mixing 執行結束
    double max = -INF;
    stFitness.reset ();

    for (int i = 0; i < nCurrent; ++i) {
        double fitness = population[i].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        stFitness.record (fitness);

    }

    if (output)
        showStatistics ();

    ++generation;
}


bool DSMGA2::shouldTerminate () {
    bool  termination = false;

    if (maxFe != -1) {
        if (Chromosome::nfe > maxFe)
            termination = true;
    }

    if (maxGen != -1) {
        if (generation > maxGen)
            termination = true;
    }

    if (population[0].getMaxFitness() <= stFitness.getMax() )
        termination = true;


    if (stFitness.getMax() - EPSILON <= stFitness.getMean() )
        termination = true;

    return termination;

}


bool DSMGA2::foundOptima () {
    return (stFitness.getMax() > population[0].getMaxFitness());
}


void DSMGA2::showStatistics () {

    printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f \n ",
            generation, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin ());
    fflush(NULL);
}



void DSMGA2::buildFastCounting() {

    if (SELECTION) {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[selectionIndex[i]].getVal(j));
            }

    } else {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[i].getVal(j));
            }
    }

}

int DSMGA2::countOne(int x) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;

        val = fastCounting[x].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countXOR(int x, int y) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;


        val = fastCounting[x].gene[i];

        val ^= fastCounting[y].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}

//2016-03-09
// Almost identical to DSMGA2::findClique
// except check 00 or 01 before adding connection
void DSMGA2::findMask(Chromosome& ch, list<int>& result,int startNode){
    result.clear();

    
	DLLA rest(ell);
	genOrderELL();
	for( int i = 0; i < ell; i++){
		if(orderELL[i] == startNode)
			result.push_back(orderELL[i]);
		else
			rest.insert(orderELL[i]);
	}

	double *connection = new double[ell];

	for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
	    pair<double, double> p = graph(startNode, *iter);
		int i = ch.getVal(startNode);
		int j = ch.getVal(*iter);
		if(i == j)//p00 or p11
			connection[*iter] = p.first;
		else      //p01 or p10
			connection[*iter] = p.second;
	}
   
    while(!rest.isEmpty()){

        // cout << "A" << endl;

	    double max = -INF;
		int index = -1;
		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
		    if(max < connection[*iter]){
                // cout << "B" << endl;
			    max = connection[*iter];
				index = *iter;
			}
		}

		rest.erase(index);
		result.push_back(index);

		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){

            if(index == -1)
            {
                break;
            }
            
            // cout << "index= " <<index<< endl;
			pair<double, double> p = graph(index, *iter);
			int i = ch.getVal(index);
			int j = ch.getVal(*iter);

			if(i == j)//p00 or p11
            {
                // cout << "C" << endl;
				connection[*iter] += p.first;
                // cout << "D" << endl;
            }
			else      //p01 or p10
            {
                // cout << "E" << endl;
				connection[*iter] += p.second;
                // cout << "F" << endl;
            }

            // if (*iter >= 0 && *iter < ell) {
            //     if (i == j) // p00 or p11
            //         connection[*iter] += p.first;
            //     else        // p01 or p10
            //         connection[*iter] += p.second;
            // } else {
            //     std::cerr << "Error: *iter out of bounds! Value: " << *iter << std::endl;
            // }
		}

        if(index == -1)
        {
            break;
        }
	}



	delete []connection;
  
}

void DSMGA2::restrictedMixing(Chromosome& ch) {
    
    int startNode = myRand.uniformInt(0, ell - 1);    
    


    list<int> mask;
	findMask(ch, mask,startNode);
    size_t size = findSize(ch, mask);
   
    list<int> mask_size; 
    findMask_size(ch,mask_size,startNode,size);
    size_t size_original = findSize(ch,mask_size);

    if (size > size_original)
        size = size_original;
    while (mask.size() > size)
        mask.pop_back();
   
    bool taken = restrictedMixing(ch, mask);


    EQ = true;
    if (taken) {
    
        genOrderN();

        for (int i=0; i<nCurrent; ++i) {
            // if (shouldBackMixing(ch, mask)) {
            if (1) {
            if (EQ)
                backMixingE(ch, mask, population[orderN[i]]);
            else
                backMixing(ch, mask, population[orderN[i]]);
            }
        }
    }

}

bool DSMGA2::shouldBackMixing(Chromosome& source, list<int>& mask) {
	bool result = true;

	for (list<int>::iterator it = mask.begin(); it != mask.end() && result; ++it) {
		int allele = source.getVal(*it);

		int differentCounter = 0;

		for (int i = 0; i < nCurrent && differentCounter < 2; i++) {
			if (population[i].getVal(*it) != allele) {
				differentCounter++;
			}
		}

		result = differentCounter >= 2 || differentCounter == 0;
	}

	return result;
}



void DSMGA2::findMask_size(Chromosome& ch, list<int>& result,int startNode,int bound){
    result.clear();

    
	DLLA rest(ell);

	for( int i = 0; i < ell; i++){
		if(orderELL[i] == startNode)
			result.push_back(orderELL[i]);
		else
			rest.insert(orderELL[i]);
	}

	double *connection = new double[ell];

	for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
	    pair<double, double> p = graph_size(startNode, *iter);
		int i = ch.getVal(startNode);
		int j = ch.getVal(*iter);
		if(i == j)//p00 or p11
			connection[*iter] = p.first;
		else      //p01 or p10
			connection[*iter] = p.second;
	}
    bound--;
    while(!rest.isEmpty()&&bound>0){
        bound--;
	    double max = -INF;
		int index = -1;
		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
		    if(max < connection[*iter]){
			    max = connection[*iter];
				index = *iter;
			}
		}

		rest.erase(index);
		result.push_back(index);

		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
			pair<double, double> p = graph_size(index, *iter);
			int i = ch.getVal(index);
			int j = ch.getVal(*iter);
			if(i == j)//p00 or p11
				connection[*iter] += p.first;
			else      //p01 or p10
				connection[*iter] += p.second;
		}
	}

    delete []connection;
  
}

void DSMGA2::backMixing(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        des = trial;
          
        return;
    }

}

void DSMGA2::backMixingE(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        EQ = false;
        des = trial;

return;
    }

    //2016-10-21
    if (trial.getFitness() >= des.getFitness() - EPSILON) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        des = trial;
        return;
    }

}

bool DSMGA2::restrictedMixing(Chromosome& ch, list<int>& mask) {

    bool taken = false;
    size_t lastUB = 0;

    for (size_t ub = 1; ub <= mask.size(); ++ub) {

        size_t size = 1;
        Chromosome trial(ell);
        trial = ch;
	    
		//2016-03-03
	    vector<int> takenMask;

        for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
            
            //2016-03-03
			takenMask.push_back(*it);

            trial.flip(*it);

            ++size;
            if (size > ub) break;
        }

        //if (isInP(trial)) continue;
        //2016-10-21
        if (isInP(trial)) break;

        if (trial.getFitness() >= ch.getFitness() - EPSILON) {
            pHash.erase(ch.getKey());
            pHash[trial.getKey()] = trial.getFitness();

            taken = true;
            ch = trial;
        }

        if (taken) {
            lastUB = ub;
            break;
        }
    }

    if (lastUB != 0) {
        while (mask.size() > lastUB)
            mask.pop_back();
    }

    return taken;

}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask) const {

    DLLA candidate(nCurrent);
    for (int i=0; i<nCurrent; ++i)
        candidate.insert(i);

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

        int allele = ch.getVal(*it);

        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            if (population[*it2].getVal(*it) == allele)
                candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;

        ++size;
    }

    return size;


}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask, Chromosome& ch2) const {

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
        if (ch.getVal(*it) == ch2.getVal(*it)) break;
        ++size;
    }
    return size;
}

void DSMGA2::mixing() {

    if (SELECTION)
        selection();

    //* really learn model
    buildFastCounting();
    buildGraph();
    buildGraph_sizecheck();
    for (int i=0; i<ell; ++i)
       findClique(i, masks[i]); // replaced by findMask in restrictedMixing

    int repeat = (ell>50)? ell/50: 1;

    for (int k=0; k<repeat; ++k) {

        genOrderN();
        for (int i=0; i<nCurrent; ++i) {
            restrictedMixing(population[orderN[i]]);
            if (Chromosome::hit) break;
        }
        if (Chromosome::hit) break;
    }


}

inline bool DSMGA2::isInP(const Chromosome& ch) const {

    unordered_map<unsigned long, double>::const_iterator it = pHash.find(ch.getKey());
    return (it != pHash.end());
}

inline void DSMGA2::genOrderN() {
    myRand.uniformArray(orderN, nCurrent, 0, nCurrent-1);
}

inline void DSMGA2::genOrderELL() {
    myRand.uniformArray(orderELL, ell, 0, ell-1);
}

void DSMGA2::buildGraph() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;

            double linkage = computeMI(p00,p01,p10,p11);
            
            //2016-04-08_computeMI_entropy
            double linkage00 = 0.0, linkage01 = 0.0;
            if (p00 > EPSILON)
                linkage00 += p00*log(p00/p_0/p0_);
            if (p11 > EPSILON)
                linkage00 += p11*log(p11/p_1/p1_);
            if (p01 > EPSILON)
                linkage01 += p01*log(p01/p0_/p_1);
            if (p10 > EPSILON)
                linkage01 += p10*log(p10/p1_/p_0);
           
            if(Chromosome::nfe < 0){
                pair<double, double> p(linkage, linkage);
                graph.write(i, j, p);
            }
            else{
                pair<double, double> p(linkage00, linkage01);
                graph.write(i, j, p);
            }
				
        }
    }


    delete []one;

}
void DSMGA2::buildGraph_sizecheck() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;

            double linkage = computeMI(p00,p01,p10,p11);
            
            //2016-04-08_computeMI_entropy
            double linkage00 = 0.0, linkage01 = 0.0;
            if (p00 > EPSILON)
                linkage00 += p00*log(p00/p_0/p0_);
            if (p11 > EPSILON)
                linkage00 += p11*log(p11/p_1/p1_);
            if (p01 > EPSILON)
                linkage01 += p01*log(p01/p0_/p_1);
            if (p10 > EPSILON)
                linkage01 += p10*log(p10/p1_/p_0);
        
	
            pair<double, double> p(linkage, linkage);
            graph_size.write(i, j, p);
			
        }
    }


    delete []one;

}


// from 1 to ell, pick by max edge
void DSMGA2::findClique(int startNode, list<int>& result) {
    
   }

    
    double DSMGA2::computeMI(double a00, double a01, double a10, double a11) const {

    double p0 = a00 + a01;
    double q0 = a00 + a10;
    double p1 = 1-p0;
    double q1 = 1-q0;

    double join = 0.0;
    if (a00 > EPSILON)
        join += a00*log(a00);
    if (a01 > EPSILON)
        join += a01*log(a01);
    if (a10 > EPSILON)
        join += a10*log(a10);
    if (a11 > EPSILON)
        join += a11*log(a11);

    double p = 0.0;
    if (p0 > EPSILON)
        p += p0*log(p0);
    if (p1 > EPSILON)
        p += p1*log(p1);


    double q = 0.0;
    if (q0 > EPSILON)
        q += q0*log(q0);
    if (q1 > EPSILON)
        q += q1*log(q1);

    return -p-q+join;

}


void DSMGA2::selection () {
    tournamentSelection ();
}


// tournamentSelection without replacement
void DSMGA2::tournamentSelection () {
    int i, j;

    int randArray[selectionPressure * nCurrent];

    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray (randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);

    for (i = 0; i < nCurrent; i++) {

        int winner = 0;
        double winnerFitness = -INF;

        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = population[challenger].getFitness();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }

        }
        selectionIndex[i] = winner;
    }
}


int DSMGA2::get_ch_dist(Chromosome x, Chromosome y) {
    
    int result(0);
    
    // for (int i=0; i!=ell; ++i) {
    //     if (x.getVal(i) != y.getVal(i))
    //         ++result;
    // }
    
    // assert(result >=0 && result <= ell);

    for (int i=0; i<ell; ++i) {
        result += myBD.countOne(x.getVal(i) ^ y.getVal(i));
    }
    
    return result;
}
