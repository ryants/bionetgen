# -*- coding: utf-8 -*-
"""
Created on Wed May 30 11:44:17 2012

@author: proto
"""
from copy import deepcopy
from lxml import etree
import pygraphviz as pgv

class Species:
    def __init__(self):
        self.molecules = []
        self.bondNumbers = []
        self.bonds = []
        
    
    
    def getBondNumbers(self):
        bondNumbers = [0]
        for element in self.molecules:
            bondNumbers.extend(element.getBondNumbers())
        return bondNumbers        
        
    def copy(self):
        species = Species()
        for molecule in self.molecules:
            species.molecules.append(molecule.copy())
        return species
    
    def addMolecule(self,molecule,concatenate=False,iteration = 1):
        if not concatenate:
            self.molecules.append(molecule)
        else:
            counter = 1
            for element in self.molecules:
                
                if element.name == molecule.name:
                    if iteration == counter:
                        element.extend(molecule)
                        return
                    else:
                        counter +=1
            self.molecules.append(molecule)
            #self.molecules.append(molecule)
            #for element in self.molecules:
            #    if element.name == molecule.name:
                
    def addCompartment(self,tags):
        for molecule in self.molecules:
            molecule.setCompartment(tags)
    
    def deleteMolecule(self,moleculeName):
        deadMolecule = None
        for element in self.molecules:
            if element.name == moleculeName:
                deadMolecule = element
        if deadMolecule == None:
            return
        bondNumbers = deadMolecule.getBondNumbers()
        self.molecules.remove(deadMolecule)
        for element in self.molecules:
            for component in element.components:
                for number in bondNumbers:
                    if str(number) in component.bonds:
                        component.bonds.remove(str(number))
    
    def getMolecule(self,moleculeName):
        for molecule in self.molecules:
            if moleculeName == molecule.name:
                return molecule
        return None
    def getSize(self):
        return len(self.molecules)
        
    def getMoleculeNames(self):
        return [x.name for x in self.molecules]
    
    def contains(self,moleculeName):
        for molecule in self.molecules:
            if moleculeName == molecule.name:
                return True
        return False
        
    def addChunk(self,tags,moleculesComponents,precursors):
        '''
        temporary transitional method
        '''
        for (tag,components) in zip (tags,moleculesComponents):
            if self.contains(tag):
                tmp = self.getMolecule(tag)
            else:
                tmp = Molecule(tag)
                #for element in precursors:
                #    if element.getMolecule(tag) != None:
                #        tmp = element.getMolecule(tag)
            
            for component in components:
                if tmp.contains(component[0][0]):
                    tmpCompo = tmp.getComponent(component[0][0])
                    #continue
                else:
                    tmpCompo = Component(component[0][0])
                
                for index in range(1,len(component[0])):
                    tmpCompo.addState(component[0][index])
                if len(component) > 1:
                    tmpCompo.addBond(component[1])
                if not tmp.contains(component[0][0]):   
                    tmp.addComponent(tmpCompo)
            if not self.contains(tag):
                self.molecules.append(tmp)
                
    def extend(self,species,update=True):
        if(len(self.molecules) == len(species.molecules)):
            for (selement,oelement) in zip(self.molecules,species.molecules):
                for component in oelement.components:
                    if component.name not in [x.name for x in selement.components]:
                        selement.components.append(component)
                    else:
                        for element in selement.components:
                            if element.name == component.name:
                                element.addStates(component.states,update)
                                
        else:
            for element in species.molecules:
                if element.name not in [x.name for x in self.molecules]:
                    
                    self.addMolecule(deepcopy(element),update)
                else:
                    for molecule in self.molecules:
                        if molecule.name == element.name:
                            for component in element.components:
                                if component.name not in [x.name for x in molecule.components]:
                                    molecule.addComponent(deepcopy(component),update)
                                else:
                                    comp = molecule.getComponent(component.name)
                                    for state in component.states:
                                        comp.addState(state,update)
                    
    
    def updateBonds(self,bondNumbers):
        newBondNumbers = deepcopy(bondNumbers)
        correspondence = {}
        intersection = [int(x) for x in newBondNumbers if x in self.getBondNumbers()]
        for element in self.molecules:
            for component in element.components:
                for index in range(0,len(component.bonds)):
                    if int(component.bonds[index]) in intersection:
                        
                        if component.bonds[index] in correspondence:
                            component.bonds[index] = correspondence[component.bonds[index]]
                        else:
                            correspondence[component.bonds[index]] = max(intersection) + 1
                            component.bonds[index] = max(intersection) + 1
                        #intersection = [int(x) for x in newBondNumbers if x in self.getBondNumbers()]
    
    def append(self,species):
        newSpecies = (deepcopy(species))
        newSpecies.updateBonds(self.getBondNumbers())
        
        for element in newSpecies.molecules:
            self.molecules.append(deepcopy(element))              
        
    def __str__(self):
        return '.'.join([x.toString() for x in self.molecules])
        
    def str2(self):
        return '.'.join([x.str2() for x in self.molecules])
        
    def reset(self):
        for element in self.molecules:
            element.reset()
            
    def toString(self):
        return self.__str__()
        
    def createGraph(self):
        speciesGraph = nx.Graph()
        for molecule in self.molecules():
            speciesGraph.add_node(molecule.createGraph())
        return speciesGraph
        
    def graphVizGraph(self,graph,identifier):
        speciesDictionary = {}

        s1 = graph.subgraph(name = "cluster%s" % identifier)
        for idx,molecule in enumerate(self.molecules):
            ident = "%s_m%i" %(identifier,idx)
            speciesDictionary[molecule.idx] = ident
            compDictionary = molecule.graphVizGraph(s1,ident)
            speciesDictionary.update(compDictionary)
        for bond in self.bonds:
            s1.add_edge(speciesDictionary[bond[0]],speciesDictionary[bond[1]],directed=False)
        return speciesDictionary
        
    def containsComponentIdx(self,idx,dictionary):
        for molecule in self.molecules:
            for component in molecule.components:
                if component.idx == idx:
                    return dictionary[idx]
        return None
        
    def notContainsComponentIdx(self,idx):
        context = []
        for molecule in self.molecules:
            for component in molecule.components:
                if component.idx not in idx:
                    context.append(component)
        return context
        

class Molecule:
    def __init__(self,name,idx):
        self.components = []
        self.name = name
        self.compartment = ''
        self.idx = idx
        
    def copy(self):
        molecule = Molecule(self.name)
        for element in self.components:
            molecule.components.append(element.copy())
        return molecule 
        
    def addChunk(self,chunk):
        component = Component(chunk[0][0][0][0])
        component.addState(chunk[0][0][0][1])
        self.addComponent(component)
        
    def addComponent(self,component,overlap=0):
        if not overlap:
            self.components.append(component)
        else:
            if not component.name in [x.name for x in self.components]:
                self.components.append(component)
            else:
                compo = self.getComponent(component.name)
                for state in component.states:
                    compo.addState(state)
    
    def setCompartment(self,compartment):
        self.compartment = compartment   
             
    def getBondNumbers(self):
        bondNumbers = []
        for element in self.components:
                bondNumbers.extend([int(x) for x in element.bonds])
        return bondNumbers
        
    def getComponent(self,componentName):
        for component in self.components:
            if componentName == component.getName():
                return component
                
    def removeComponent(self,componentName):
        x = [x for x in self.components if x.name == componentName]
        if x != []:
            self.components.remove(x[0])
            
    def removeComponents(self,components):
        for element in components:
            if element in self.components:
                self.components.remove(element)
                
    def addBond(self,componentName,bondName):
        bondNumbers = self.getBondNumbers()
        while bondName in bondNumbers:
            bondName += 1
        component = self.getComponent(componentName)
        component.addBond(bondName)
        
    def getComponentWithBonds(self):
        return [x for x in self.components if x.bonds != []]
        
    def contains(self,componentName):
        return componentName in [x.name for x in self.components]
        
    def __str__(self):
        self.components.sort()
        return self.name + '(' + ','.join([str(x) for x in self.components]) + ')' + self.compartment
        
    def toString(self):
        return self.__str__()
        
    def str2(self):
        self.components.sort()
        return self.name + '(' + ','.join([x.str2() for x in self.components]) + ')'
        
    def extend(self,molecule):
        for element in molecule.components:
            comp = [x for x in self.components if x.name == element.name]
            if len(comp) == 0:
                self.components.append(deepcopy(element))
            else:
                for bond in element.bonds:
                    comp[0].addBond(bond)
                for state in element.states:
                    comp[0].addState(state)
                    
    def reset(self):
        for element in self.components:
            element.reset()
            
    def update(self,molecule):
        for comp in molecule.components:
            if comp.name not in [x.name for x in self.components]:
                self.components.append(deepcopy(comp))
                
    def graphVizGraph(self,graph,identifier,components=None):
        moleculeDictionary = {}
        if len(self.components) == 0:
            graph.add_node(identifier,label=self.name)
            moleculeDictionary[self.idx] = identifier
        else:
            s1 = graph.subgraph(name = "cluster%s" % identifier,label=self.name)
            if components == None:
                tmpComponents = self.components
            else:
                tmpComponents = components
            for idx,component in enumerate(tmpComponents):
                ident = "%s_c%i" %(identifier,idx)
                moleculeDictionary[component.idx] = ident
                compDictionary = component.graphVizGraph(s1,ident)
                moleculeDictionary.update(compDictionary)
        return moleculeDictionary
                
    def createGraph(self,identifier):
        molDictionary = {}
        molecule = etree.Element("graph")
        molecule.set("id","m%s:" % identifier)
        molecule.set("name",self.name)
        for idx,component in enumerate(self.components):
            gcomponent = etree.Element("node")
            gcomponent.set("id","m%s::c%s" % (identifier,idx))
            print component.createGraph("m%s::c%s" % (identifier,idx))
            compGraph,compDictionary = component.createGraph("m%s::c%s" % (identifier,idx))
            molDictionary.update(compDictionary)
            molDictionary[component.idx] = gcomponent.get("id")
            if compGraph != None:
                gcomponent.append(compGraph)
            else:
                gcomponent.set("name",component.name)
            molecule.append(gcomponent)
        
        '''
        moleculeGraph = nx.Graph(name = self.name)
        for component in self.components:
            moleculeGraph.add_node(component.createGraph())
        return moleculeGraph
        '''
        return molecule,molDictionary       
    
class Component:
    def __init__(self,name,idx,bonds = [],states=[]):
        self.name = name
        self.states = []
        self.bonds = []
        self.activeState = ''
        self.idx = idx
        
    def copy(self):
        component = Component(self.name,deepcopy(self.bonds),deepcopy(self.states))
        component.activeState = deepcopy(self.activeState)     
        return component
        
    def addState(self,state,update=True):
        if not state in self.states:
            self.states.append(state)
        if update:
            self.setActiveState(state)
        #print 'LALALA',state
    def addStates(self,states,update=True):
        for state in states:
            if state not in self.states:
                self.addState(state,update)
        
    def addBond(self,bondName):
        if not bondName in self.bonds:
            self.bonds.append(bondName)
        
    def setActiveState(self,state):
        if state not in self.states:
            return False
        self.activeState = state
        return True
        
    def getRuleStr(self):
        tmp = self.name
        if len(self.bonds) > 0:
            tmp += '!' + '!'.join([str(x) for x in self.bonds])
        if self.activeState != '':
            tmp += '~' + self.activeState
        return tmp
        
    def getTotalStr(self):
        return self.name + '~'.join(self.states)
    
    def getName(self):
        return self.name 
        
    def __str__(self):
        return self.getRuleStr()
        
    def str2(self):
        tmp = self.name
        if len(self.bonds) > 0:
            tmp += '!' + '!'.join([str(x) for x in self.bonds])
        if len(self.states) > 0:
            tmp += '~' + '~'.join([str(x) for x in self.states])
        return tmp        
        
    def __hash__(self):
        return self.name
        
    def reset(self):
        self.bonds = []
        if 'U' in self.states:
            self.activeState = 'U'
    
    def graphVizGraph(self,graph,identifier):
        compDictionary = {}
        if len(self.states) == 0:
            graph.add_node(identifier,label=self.name)
        else:
            s1 = graph.subgraph(name="cluster%s" % identifier,label=self.name)
            if self.activeState != '':
                ident = "%s_s%i" % (identifier,0)
                compDictionary[self.activeState] = ident
                s1.add_node(ident,label = self.activeState)
            else:
                for idx,state in enumerate(self.state):
                    ident = "%s_s%i" % (identifier,idx)
                    s1.add_node(ident,label = state)
                    compDictionary[state] = ident
        return compDictionary
            
    def createGraph(self,identifier):
        '''
        componentGraph = nx.Graph(name = 'self.name')
        for state in self.states:        
            componentGraph.add_node(state)
        '''
        compDictionary = {}
        if len(self.states) == 0:
            return None,{}
        component = etree.Element("graph")
        component.set("id","%s:" % identifier)
        component.set("name",self.name)
        for idx,element in enumerate(self.states):
            state = etree.Element("node")
            state.set("id",("%s::s%i") % (identifier,idx))
            state.set("name",element)
            #TODO: adjust so that you can handle simmetry
            compDictionary[state.get("id")] = element 
            component.append(state)
        
        
        return component,compDictionary
        
class Action:
    def __init__(self):
        self.action = ''
        self.site1 = ''
        self.site2 = ''
    
    def setAction(self,action,site1,site2=''):
        self.action = action
        self.site1 = site1
        self.site2 = site2
        
    def __str__(self):
        return '%s, %s, %s' %(self.action,self.site1,self.site2)
        
class Databases:
    def __init__(self):
        self.translator ={}
        self.synthesisDatabase = {}
        self.catalysisDatabase = {}
        self.rawDatabase = {}
        self.labelDictionary = {}
        self.synthesisDatabase2 = {}
        
    def getRawDatabase(self):
        return self.rawDatabase
        
    def getLabelDictionary(self):
        return self.labelDictionary
        
    def add2LabelDictionary(self,key,value):
        temp = tuple(key)
        temp = temp.sort()
        self.labelDictionary[temp] = value

    def add2RawDatabase(self,rawDatabase):
        pass
    
    def getTranslator(self):
        return self.translator
    