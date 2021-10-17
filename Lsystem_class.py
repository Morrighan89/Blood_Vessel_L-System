"""
MIT License

Copyright (c) 2021 Riccardo Ferrero Scott Chacon and others 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


import random
import math
import numpy as np
from tqdm import tqdm

class Lsystem:
    def __init__(self, axiom=[], ruleset={},alphabet=''):
        self.axiom = axiom
        self.ruleset = ruleset
        self.generation = 0
        self.alphabet=alphabet

    def process(self, original_string):
        #  process a string with ruleset and return next gen string
        tranformed_string = []
        for letter in original_string:
            if letter[0] in self.alphabet:
                currentNode=node(letter[0],letter[1],letter[2])
                tranformed_string.extend(ruleOutput(currentNode,self.ruleset[letter[0]]))
            else:
                tranformed_string.append(letter)
        self.generation += 1
        return tranformed_string

    def processGen(self, gen=1):
        #  process strings until the given generation
        actualString = self.axiom
        for i in tqdm(range(gen)):
            actualString = self.process(actualString)
        return actualString

    def getGeneration(self):
        return self.generation

class rule:
    def __init__(self,ruleType,proba=1):
        self.ruleTypes=['fw','fwturn','bif','end']
        self.proba = proba
        self.ruleType=ruleType
        if ruleType not in self.ruleTypes:
            print("Type %s does not exist, please check it" % ruleType)
            return

    def evaluate(self,node):
        rulevalues=[]
        if self.ruleType == 'fw':
            rulevalues.append(['f',node.diam,node.len])
            rulevalues.append([node.vessKind,node.diam,node.len])
            return rulevalues
        elif self.ruleType == 'bif':
            return writeBifurcation(node)
        elif self.ruleType == 'fwturn':
            if random.random()<0.5:
                rulevalues.append(['*',random.uniform(0,np.pi/4)])
            else:
                rulevalues.append(['/',random.uniform(0,np.pi/4)])
            if random.random()<0.5:
                rulevalues.append(['+',random.uniform(0,np.pi/4)])
            else:
                rulevalues.append(['-',random.uniform(0,np.pi/4)])
            rulevalues.append(['f',node.diam,node.len])
            rulevalues.append([node.vessKind,node.diam,node.len])
        elif self.ruleType == 'end':
            rulevalues.append(['E',0,0])
        return rulevalues


def checkRules(rules=[]):
    sum = 0
    for r in rules:
        sum += r.proba
    return True if sum == 1 else False

def ruleOutput(node,rules=[]):
    if checkRules(rules):
        rand_num = random.uniform(0,1)
        sum = 0
        for r in rules:
            sum += r.proba
            if rand_num < sum:
                return r.evaluate(node)
    else:
        raise Exception('the sum of probabilities is expected to be 1')

class node:
    def __init__(self,vessKind='A',diam=1,len=1):
        self.vessKind=vessKind
        self.diam = diam
        self.len = len
    def getNodeParam(self):
        return '('+str(self.diam)+','+str(self.len)+')'
def writeBifurcation(node):
    ruleString = []
    turning=random.random()
    angle=random.uniform(0,np.pi/2)
    ruleString.append(['f',node.diam,node.len])
    params=calculateBifurcation(node.vessKind,node.diam,node.len,random.uniform(0.1,1))
    ruleString.append('[')
    if turning<0.3:
        ruleString.append(['*',angle])
    elif turning<0.6:
        ruleString.append(['/',angle])
    else:
        None
    ruleString.append(['+',params['th1']])
    node.vessKind=vesselKindEvaluation(params['d1'])
    ruleString.append([node.vessKind,params['d1'],params['l1']])
    ruleString.append(']')
    ruleString.append('[')
    if turning<0.3:
        ruleString.append(['*',angle])
    elif turning<0.6:
        ruleString.append(['/',angle])
    else:
        None
    ruleString.append(['-',params['th2']])
    node.vessKind=vesselKindEvaluation(params['d2'])
    ruleString.append([node.vessKind,params['d2'],params['l2']])
    ruleString.append(']')
    return ruleString

def calculateBifurcation(kind,d0,l0,alpha=1):
    params=dict.fromkeys(('d1','l1','th1','ph1','d2','l2','th2','ph2'))
    coin=random.random()
    if kind == 'A':
        alpha=alpha*0.7
        gamma1=1/(1+alpha**3)**(1/3)
        gamma2=alpha/(1+alpha**3)**(1/3)
        lambda1=gamma1
        lambda2=gamma2
    elif kind == 'B':
        alpha=alpha*0.8
        gamma1=1/(1+alpha**3)**(1/3)
        gamma2=alpha/(1+alpha**3)**(1/3)
        lambda1=gamma1
        lambda2=gamma2
    else:
        gamma1=1/(1+alpha**3)**(1/3)
        gamma2=alpha/(1+alpha**3)**(1/3)
        lambda1=1
        lambda2=1
    th1=((1+alpha**3)**(4/3)+1-alpha**4)/(2*(1+alpha**3)**(2/3))
    th2=((1+alpha**3)**(4/3)+alpha**4-1)/(2*alpha**2*(1+alpha**3)**(2/3))
    if coin < 0.5:
        params['d1']=d0*gamma1
        params['d2']=d0*gamma2
        params['l1']=l0*lambda1
        params['l2']=l0*lambda2
        params['th1']=math.acos(th1)
        params['th2']=math.acos(th2)
    else:
        params['d1']=d0*gamma2
        params['d2']=d0*gamma1
        params['l1']=l0*lambda2
        params['l2']=l0*lambda1
        params['th1']=math.acos(th2)
        params['th2']=math.acos(th1)
    params['ph1']=random.random()
    params['ph2']=params['ph1']
    return params

def vesselKindEvaluation(diam):
    if diam<1/100:
        vessKind='E'
    elif diam<1/50:
        vessKind='C'
    elif diam<1/4:
        vessKind='B'
    else:
        vessKind='A'
    return vessKind
