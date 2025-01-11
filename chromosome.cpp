/***************************************************************************
 *   Copyright (C) 2015 by TEIL                                            *
 ***************************************************************************/

#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include "spin.h"
#include "chromosome.h"
#include "nk-wa.h"
#include "sat.h"
#include "dsmga2.h"
#include <iostream>
#define TRAP_K 4
#define TRAP_K_3 3

Chromosome::Chromosome()
{
    length = 0;
    lengthLong = 0;
    gene = NULL;
    evaluated = false;
}

Chromosome::Chromosome(int n_length)
{
    gene = NULL;
    init(n_length);
}

Chromosome::~Chromosome()
{
    if (gene != NULL)
        delete[] gene;
}

void Chromosome::init(int _length)
{
    length = _length;
    lengthLong = quotientLong(length) + 1;

    if (gene != NULL)
        delete[] gene;

    gene = new unsigned long[lengthLong];
    gene[lengthLong - 1] = 0;

    evaluated = false;
}

void Chromosome::init0(int _length)
{
    length = _length;
    lengthLong = quotientLong(length) + 1;

    if (gene != NULL)
        delete[] gene;

    gene = new unsigned long[lengthLong];

    for (int i = 0; i < lengthLong; ++i)
        gene[i] = 0;

    key = 0;
    evaluated = false;
}

void Chromosome::initR(int _length)
{
    length = _length;
    lengthLong = quotientLong(length) + 1;

    if (gene != NULL)
        delete[] gene;

    gene = new unsigned long[lengthLong];
    gene[lengthLong - 1] = 0;

    key = 0;
    for (int i = 0; i < length; ++i)
    {

        int val = myRand.flip();
        setValF(i, val);
        if (val == 1)
            key ^= zKey[i];
    }

    evaluated = false;
}

double Chromosome::getFitness()
{   
    // cout << "hihihi+++++++++++++++++++" << endl;
    if (evaluated){
        // cout << "evaluated+++++++++++++" << endl;
        return fitness;
    }
    else
    {
        // cout << "not evaluated+++++++++++++" << endl;
        fitness = evaluate();
        if (!hit && fitness > getMaxFitness())
        {

            // cout<< "Generation: " << generation << " ";
            cout << "最佳解!!!!!" << endl;

            for (int i = 0; i < length; ++i)
                cout << getVal(i);
            cout << endl;

            hit = true;
            hitnfe = nfe + lsnfe;
        }
        return fitness;
    }
}

bool Chromosome::isEvaluated() const
{
    return evaluated;
}

bool Chromosome::hasSeen() const
{

    unordered_map<unsigned long, double>::iterator it = cache.find(key);
    if (it != cache.end())
        return true;
    return false;
}

double Chromosome::evaluate()
{

    if (CACHE)
        if (hasSeen())
        {
            evaluated = true;
            return cache[key];
        }

    ++nfe;
    evaluated = true;
    double accum = 0.0;

    switch (function)
    {
    case ONEMAX:
        accum = oneMax();
        break;
    case MKTRAP:
        accum = mkTrap(1, 0.8);
        break;
    case CYCTRAP:
        accum = cycTrap(1, 0.8);
        break;
    case FTRAP:
        accum = fTrap();
        break;
    case SPINGLASS:
        accum = spinGlass();
        break;
    case NK:
        accum = nkFitness();
        break;
    case SAT:
        accum = satFitness();
        break;
    case mixTRAP:
        accum = mixTrap(4.0, 3.0, 3);
        break;
    case THREEOPTIMUM:
        accum = threeOptimum(3.0, 2.0);
        break;
    default:
        accum = mkTrap(1, 0.8);
        break;
    }

    if (CACHE)
        cache[key] = accum;

    return accum;
}

double
Chromosome::spinGlass() const
{

    int *x = new int[length];
    double result;

    for (int i = 0; i < length; i++)
        if (getVal(i) == 1)
            x[i] = 1;
        else
            x[i] = -1;

    result = evaluateSPIN(x, &mySpinGlassParams);

    delete[] x;

    return result;
}

double Chromosome::nkFitness() const
{
    char *x = new char[length];

    for (int i = 0; i < length; ++i)
    {
        x[i] = (char)getVal(i);
    }

    double result = evaluateNKProblem(x, &nkwa);
    // double result = evaluateNKWAProblem(x, &nkwa);
    delete[] x;
    return result;
}

// OneMax
double Chromosome::oneMax() const
{

    double result = 0;

    for (int i = 0; i < length; ++i)
        result += getVal(i);

    return result;
}

bool Chromosome::operator==(const Chromosome &c) const
{
    if (length != c.length)
        return false;

    for (int i = 0; i < lengthLong; i++)
        if (gene[i] != c.gene[i])
            return false;

    return true;
}

Chromosome &Chromosome::operator=(const Chromosome &c)
{

    if (length != c.length)
    {
        length = c.length;
        init(length);
    }

    evaluated = c.evaluated;
    fitness = c.fitness;
    lengthLong = c.lengthLong;
    key = c.key;

    memcpy(gene, c.gene, sizeof(long) * lengthLong);

    return *this;
}

double Chromosome::trap(int unitary, double fHigh, double fLow, int trapK) const
{
    if (unitary > trapK)
        return 0;

    if (unitary == trapK)
        return fHigh;
    else
        return fLow - unitary * fLow / (trapK - 1);
}

double Chromosome::trap_zeromax(int unitary, double fHigh, double fLow, int trapK) const
{
    if (unitary > trapK)
        outputErrMsg("unitary > trapK");
    if (unitary == 0)
        return fHigh;
    if (unitary == trapK)
        return fLow;
    else
        return fLow - ((trapK - unitary) * fLow / (trapK - 1));
}

double Chromosome::fTrap() const
{

    double result = 0.0;

    for (int i = 0; i < length / 6; ++i)
    {
        int u = 0;
        for (int j = 0; j < 6; ++j)
            u += getVal(i * 6 + j);

        if (u == 0)
            result += 1.0;
        else if (u == 1)
            result += 0.0;
        else if (u == 2)
            result += 0.4;
        else if (u == 3)
            result += 0.8;
        else if (u == 4)
            result += 0.4;
        else if (u == 5)
            result += 0.0;
        else // u == 6
            result += 1.0;
    }

    return result;
}

double Chromosome::cycTrap(double fHigh, double fLow) const
{
    int i, j;
    int u;
    int TRAP_M = length / (TRAP_K - 1);
    if (length % (TRAP_K - 1) != 0)
        outputErrMsg("TRAP_k doesn't divide length for Cyclic Setting");
    double result = 0;
    for (i = 0; i < TRAP_M; i++)
    {
        u = 0;
        int idx = i * TRAP_K - i;
        for (j = 0; j < TRAP_K; j++)
        {
            int pos = idx + j;
            if (pos == length)
                pos = 0;
            else if (pos > length)
                outputErrMsg("CYCLIC BUG");
            //
            u += getVal(pos);
        }
        result += trap(u, fHigh, fLow, TRAP_K);
    }
    return result;
}

double Chromosome::mkTrap(double fHigh, double fLow) const
{
    int i, j;
    int u;

    int TRAP_M = length / TRAP_K;

    if (length % TRAP_K != 0)
        outputErrMsg("TRAP_K doesn't divide length");

    double result = 0;

    for (i = 0; i < TRAP_M; i++)
    {
        u = 0;
        for (j = 0; j < TRAP_K; j++)
            u += getVal(i * TRAP_K + j);

        result += trap(u, fHigh, fLow, TRAP_K);
    }

    return result;
}

int Chromosome::getLength() const
{
    return length;
}

double Chromosome::getMaxFitness() const
{

    double maxF;

    switch (function)
    {
    case ONEMAX:
        maxF = length;
        break;
    case MKTRAP:
        maxF = length / TRAP_K;
        break;
    case FTRAP:
        maxF = length / 6;
        break;
    case CYCTRAP:
        maxF = length / (TRAP_K - 1);
        break;
    case SPINGLASS:
        maxF = mySpinGlassParams.opt;
        break;
    case NK:
        maxF = nkwa.maxF;
        break;
    case SAT:
        maxF = 0;
        break;
    case mixTRAP:
        maxF = length;
        break;
    case THREEOPTIMUM:
        maxF = length;
        break;
    default:
        // Never converge
        maxF = INF;
    }

    return maxF - EPSILON;
}

// contribute to lsnfe
bool Chromosome::tryFlipping(int index)
{

    int oldNFE = nfe;

    double oldF = getFitness();
    flip(index);

    // 2016-10-21
    if (getFitness() - EPSILON <= oldF)
    {
        // if (getFitness() <= oldF) {
        flip(index);
        evaluated = true;
        fitness = oldF;

        lsnfe += nfe - oldNFE;
        nfe = oldNFE;

        return false;
    }
    else
    {

        lsnfe += nfe - oldNFE;
        nfe = oldNFE;

        return true;
    }
}

bool Chromosome::GHC()
{

    int *order = new int[length];
    myRand.uniformArray(order, length, 0, length - 1);

    bool flag = false;
    for (int i = 0; i < length; ++i)
    {
        if (tryFlipping(order[i]))
            flag = true;
    }

    delete[] order;
    return flag;
}

double Chromosome::satFitness() const
{
    int *x = new int[length];

    for (int i = 0; i < length; ++i)
    {
        x[i] = getVal(i);
    }

    double result = evaluateSAT(x, &mySAT);
    delete[] x;
    return result;
}

double Chromosome::mixTrap(double fHigh, double fLow, int overlap) const
{
    int TRAP_M = length / TRAP_K;

    if (length % TRAP_K != 0)
        outputErrMsg("TRAP_K doesn't divide length");

    if (overlap < 0 || overlap >= TRAP_K)
        outputErrMsg("Invalid overlap value");
    double result_trap1 = 0.0;
    double result_trap2 = 0.0;
    for (int i = 0; i < TRAP_M; i++)
    {
        int u1 = 0;
        for (int j = 0; j < TRAP_K; j++)
        {
            // Get the value of the bit in the current block
            u1 += getVal(i * TRAP_K + j);
        }
        result_trap1 += trap(u1, fHigh, fLow, TRAP_K);

        int u2 = 0;
        for (int j = 0; j < TRAP_K; j++)
        {
            // Adjust the starting index of the block based on overlap
            int index = (i * TRAP_K + j + overlap) % length;
            // Get the value of the bit at the adjusted index
            u2 += getVal(index);
        }
        // Add the zeromax fitness for this block
        result_trap2 += trap_zeromax(u2, fHigh, fLow, TRAP_K);
    }
    double result = std::max(result_trap1, result_trap2);
    return result;
}

double Chromosome::threeOptimum(double fHigh, double fLow) const
{
    // cout << "ThreeOptimum" << endl;
    int TRAP_M = length / TRAP_K_3;

    if (length % TRAP_K_3 != 0){
        outputErrMsg("Here");
        outputErrMsg("TRAP_K doesn't divide length!");
    }
    double result_opt1 = 0.0; // Fitness for optimal solution "100"
    double result_opt2 = 0.0; // Fitness for optimal solution "010"
    double result_opt3 = 0.0; // Fitness for optimal solution "001"

    for (int i = 0; i < TRAP_M; i++)
    {
        int blockStart = i * TRAP_K_3;

        // Calculate fitness for each of the three optimal solutions
        result_opt1 += trap_threeHamming(blockStart, "100", fHigh, fLow);
        result_opt2 += trap_threeHamming(blockStart, "010", fHigh, fLow);
        result_opt3 += trap_threeHamming(blockStart, "001", fHigh, fLow);
    }

    // Return the maximum fitness among the three optimal solutions
    return std::max(result_opt1, std::max(result_opt2, result_opt3));
}

double Chromosome::trap_threeHamming(int blockStart, const std::string &optimal, double fHigh, double fLow) const
{
    // Extract the 3-bit block as a string
    std::string block = "";
    for (int i = 0; i < TRAP_K_3; i++)
    {
        block += (getVal(blockStart + i) ? '1' : '0');
    }

    // Calculate the Hamming distance between the block and the optimal solution
    double hammingDistance = 0.0;
    for (int i = 0; i < TRAP_K_3; i++)
    {
        if (block[i] != optimal[i])
            hammingDistance++;
    }

    // Assign fitness based on Hamming distance
    if (hammingDistance == 0)
    {
        return 3.0;
    }
    else
    {
        return hammingDistance - 1.0;
    }
}

/*

*/