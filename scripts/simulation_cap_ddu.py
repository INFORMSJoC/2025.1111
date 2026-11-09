# -*- coding: utf-8 -*-
import numpy as np
import gurobipy as gu
import math
import timeit
import math
import multiprocessing
import argparse
from csv import writer

BUDGET = 100
TIME_LIMIT = 18000
ScenNum = 1000

parser = argparse.ArgumentParser(description='Job number')
parser.add_argument('--job', '-d', metavar='the-no-of-data', dest='job_num', type=int, help='Choose job number', default=0)
args = parser.parse_args()

JOB=args.job_num

JOBlist=[]
for r in [100000,200000,300000,400000]:
    JOBlist.append([r])
            
                
lambda_factor   = 0.25
DDU             = True
CAP_DDU         = True
XR_R            = JOBlist[JOB][0]
lambda_rate     = 200
capital         = 50000000
profit          = 3.8*1.2
function        = 25
DISTRIBUTION    = "Truncated_Normal"
BUDGET          = 100
deviation       = 1
kappa           = 0
sup_range       = [0.75, 1.25]
C_IMPORT        = 1.3
epsilon_factor  = 0.10
a_skew          = -1

np.random.seed(0)

num_cores =multiprocessing.cpu_count()
Distance_Dict={(0, 6): 22.7482, (0, 7): 6.1988, (0, 8): 47.561800000000005, (0, 9): 90.7128, (0, 10): 95.52260000000001, (0, 11): 97.36489999999999, (0, 12): 100.479, (0, 13): 37.6286, (0, 14): 21.967299999999998, (0, 15): 89.8982, (0, 16): 100.83489999999999, (0, 17): 48.5474, (0, 18): 22.2264, (1, 6): 0.7773, (1, 7): 19.7411, (1, 8): 33.064099999999996, (1, 9): 76.2151, (1, 10): 101.2945, (1, 11): 91.82469999999999, (1, 12): 94.9388, (1, 13): 39.942800000000005, (1, 14): 0.0036, (1, 15): 84.358, (1, 16): 106.6068, (1, 17): 45.4215, (1, 18): 0.2555, (2, 6): 80.4085, (2, 7): 92.5659, (2, 8): 50.7131, (2, 9): 4.149100000000001, (2, 10): 100.80080000000001, (2, 11): 34.4789, (2, 12): 39.3302, (2, 13): 66.3565, (2, 14): 79.62769999999999, (2, 15): 43.5084, (2, 16): 106.1131, (2, 17): 56.402, (2, 18): 79.88669999999999, (3, 6): 40.3404, (3, 7): 39.2141, (3, 8): 31.2133, (3, 9): 59.0013, (3, 10): 63.0197, (3, 11): 62.824400000000004, (3, 12): 65.9385, (3, 13): 2.4993000000000003, (3, 14): 39.5595, (3, 15): 55.357699999999994, (3, 16): 68.332, (3, 17): 13.6, (3, 18): 39.818599999999996, (4, 6): 35.6526, (4, 7): 47.81, (4, 8): 4.4589, (4, 9): 47.974599999999995, (4, 10): 96.1724, (4, 11): 65.3661, (4, 12): 68.4802, (4, 13): 34.820699999999995, (4, 14): 34.8718, (4, 15): 57.8994, (4, 16): 101.4847, (4, 17): 37.994099999999996, (4, 18): 35.1308, (5, 6): 25.1516, (5, 7): 5.206300000000001, (5, 8): 49.965300000000006, (5, 9): 93.11619999999999, (5, 10): 97.7695, (5, 11): 99.6118, (5, 12): 102.7258, (5, 13): 39.8755, (5, 14): 24.3708, (5, 15): 92.1451, (5, 16): 103.0818, (5, 17): 50.7943, (5, 18): 24.6298}

class Node:
    def __init__(self, node_id, x, y, role,demand):
        self.node_id = node_id
        self.Coord_x = x
        self.Coord_y = y
        self.role = role
        self.outgoing =[]
        self.incoming=[]
        self.cap=1000000    
        self.hd={}
        self.hprop={}
        #initial capacity of the supply node at the beginning
        self.In_cap=0
        self.w_min=0
        self.w_max=float('inf')
        self.w_0=0
        
        self.P_vec={}
        
        self.beta_max=500000
        self.Support={}
        self.mu={}
        self.Lambda={}
        self.epsilon={}
        self.M=10000000

        
        self.h={} ; self.nu={} ; self.tau={} ; self.CS_S1={} ;  self.CS_S2={} ;  self.D_S1={} ;  self.D_SD={}
        
        self.v={} ; self.ups={} ; self.CS_I1={} ; self.D_I1= {} ; self.D_ID={}

        self.xi={} ; self.psi={} ; self.CS_D1={}; self.w={}
        self.lowerR={}
        self.upperR={}

        self.Demand=demand
         
    
    def AddSupplyVariables_SP(self,model,t):
        self.h[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="h"+str(self.node_id)+","+str(t) )    
        self.nu[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="nu"+str(self.node_id)+","+str(t))    
        self.tau[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="tau"+str(self.node_id)+","+str(t) )    
    
        self.CS_S1[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="CS_S1"+str(self.node_id)+","+str(t) )    
        self.CS_S2[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="CS_S2"+str(self.node_id)+","+str(t) )    

        self.D_S1[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="D_S1"+str(self.node_id)+","+str(t) )            
        self.D_SD[t]=model.addVars(demand_count,lb=0,vtype=gu.GRB.BINARY,name=["D_SD"+str(self.node_id)+","+str(t)+","+str(k) for k in range(demand_count) ])    


    def AddSupplyVariables_MP(self,model, fixed = False):
        if not fixed:
            self.x=model.addVars(T,vtype=gu.GRB.BINARY,name=["x"+str(self.node_id)+","+str(t) for t in range(T)])
            self.y=model.addVars(T,lb=0,ub=self.cap,vtype=gu.GRB.CONTINUOUS,name=["y"+str(self.node_id)+","+str(t) for t in range(T)])
        if CAP_DDU:
            self.e=model.addVars(T,R,vtype=gu.GRB.BINARY,name=["e"+","+str(self.node_id)+","+str(t)+","+str(r) for t in range(T) for r in range(R) ])

        self.MP_h={}
    
    def AddSupplyVariables_Dual(self,model):
        self.nu=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["nu"+str(self.node_id)+","+str(t) for t in range(T)])    
        self.tau=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["tau"+str(self.node_id)+","+str(t) for t in range(T)])    


    def AddSupplyVariables_Deterministic(self,model):
        self.x=model.addVars(T,vtype=gu.GRB.BINARY,name=["x"+str(self.node_id)+","+str(t) for t in range(T)])
        self.y=model.addVars(T,lb=0,ub=self.cap,vtype=gu.GRB.CONTINUOUS,name=["y"+str(self.node_id)+","+str(t) for t in range(T)])
        
        self.h=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["h"+str(self.node_id)+","+str(t) for t in range(T)] )
    
    def AddSupplyVariables_Simulation(self,model, s):        
        self.h[s]=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["h"+str(self.node_id)+","+str(t)+f",{s}" for t in range(T)])

              
    def AddPortVariables_SP(self,model,t):
        self.v[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="v"+str(self.node_id)+","+str(t)  )
        self.ups[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="ups"+str(self.node_id)+","+str(t)  )
        
        self.CS_I1[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="CS_I1"+str(self.node_id)+","+str(t) )    
        self.D_I1[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="D_I1"+str(self.node_id)+","+str(t) )    
        self.D_ID[t]=model.addVars(demand_count,lb=0,vtype=gu.GRB.BINARY,name=["D_I1"+str(self.node_id)+","+str(t)+","+str(k) for k in range(demand_count) ]) 

    def AddPortVariables_MP(self,model):
        self.MP_v={}
        
    def AddPortVariables_Dual(self,model):
        self.ups=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["ups"+str(self.node_id)+","+str(t)  for t in range(T)])
        
    def AddPortVariables_Deterministic(self,model):
        self.v=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["v"+str(self.node_id)+","+str(t) for t in range(T)])
    
    def AddPortVariables_Simulation(self,model, s):
        self.v[s]=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["v"+str(self.node_id)+","+str(t)+f",{s}" for t in range(T)])

                
    def AddDemandVariables_MP(self,model, fixed = False):
        self.beta1=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["beta1"+","+str(self.node_id)+","+str(t) for t in range(T) ])
        self.beta2=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["beta2"+","+str(self.node_id)+","+str(t) for t in range(T) ])
        self.phi1_variables={}
        self.phi2_variables={}
        if not fixed and not CAP_DDU:
            for line in self.incoming:
                if line.From in supply_nodes:
                    self.phi1_variables[line.From] = model.addVars(T,lb=0, vtype=gu.GRB.CONTINUOUS, name= ["phi1"+","+str(self.node_id)+","+str(t) for t in range(T) ])
                    self.phi2_variables[line.From] = model.addVars(T,lb=0, vtype=gu.GRB.CONTINUOUS,  name= ["phi2"+","+str(self.node_id)+","+str(t) for t in range(T) ])
        if CAP_DDU:
            self.zeta1_variables={}
            self.zeta2_variables={}
            for line in self.incoming:
                if line.From in supply_nodes:
                    self.zeta1_variables[line.From] = model.addVars(T,R,lb=0, vtype=gu.GRB.CONTINUOUS, name= ["zeta1"+","+str(self.node_id)+","+str(t) for r in range(R)  for t in range(T) ])
                    self.zeta2_variables[line.From] = model.addVars(T,R,lb=0, vtype=gu.GRB.CONTINUOUS,  name= ["zeta2"+","+str(self.node_id)+","+str(t) for r in range(R) for t in range(T) ])


    def AddDemandVariables_SP(self,model,t):
        self.xi[t]=model.addVar(lb=self.Support[t][0],ub=self.Support[t][-1], vtype=gu.GRB.CONTINUOUS, name="xi"+","+str(self.node_id)+","+str(t))
        self.psi[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="psi"+str(self.node_id)+","+str(t) )

        self.CS_D1[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="CS_D1"+str(self.node_id)+","+str(t))    
        
    def AddDemandVariables_Dual(self,model):
        self.xi={}
        for t in range(T):
            self.xi[t]=model.addVar(lb=self.Support[t][0],ub=self.Support[t][-1], vtype=gu.GRB.CONTINUOUS,name="xi"+","+str(self.node_id)+","+str(t))

        self.psi=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["psi"+str(self.node_id)+","+str(t)  for t in range(T)])
        
        
    def AddSupplyConstraints_SP(self,model,t):
        #11f

        model.addConstr(-self.h[t]>=0 , name="P_Cap"+str(self.node_id)+','+str(t))
        #11e supply
        model.addConstr(gu.quicksum(line.z[t] for line in self.outgoing) == self.h[t])
       
        #dual constraint 1
        model.addConstr((-self.nu[t]+self.tau[t] <= self.c_VP[t] ), name='Dual_1')
        #dual constraints 3
        model.addConstrs((-self.tau[t]+self.outgoing[k].To.psi[t] <= self.outgoing[k].c_VT[t]*self.outgoing[k].dist for k in range(demand_count) ), name='Dual_3' )
        
        #complemantary slackness1

        model.addConstr((-self.h[t] <= self.M*self.CS_S1[t]),name='CS_S1C'+str(self.node_id)+','+str(t))
        model.addConstr(self.nu[t] <= self.M*(1-self.CS_S1[t]) )
        #complemantary slackness2
        model.addConstr((self.h[t]-gu.quicksum(line.z[t] for line in self.outgoing) <= self.M*self.CS_S2[t] ),name='CS_2' )
        model.addConstr((self.tau[t]<=self.M*(1-self.CS_S2[t]) ),name='CS_2')

        #dual complementart slackness 1
        model.addConstr(self.c_VP[t]+self.nu[t]-self.tau[t]<=self.M*self.D_S1[t] )
        model.addConstr(self.h[t]<=self.M*(1-self.D_S1[t]) )
    
        #dual complementart slackness 3
        model.addConstrs(self.outgoing[k].c_VT[t]*self.outgoing[k].dist+self.tau[t]-self.outgoing[k].To.psi[t]<=self.M*self.D_SD[t][k] for k in range(demand_count)  )
        model.addConstrs(self.outgoing[k].z[t]<=self.M*(1-self.D_SD[t][k]) for k in range(demand_count) )

    def AddSupplyConstraints_MP(self,model):
        #11b
        model.addConstrs(self.y[t]<=self.cap*self.x[t]  for t in range(T))
        #11c
        model.addConstrs(self.x[t]<=self.x[t+1] for t in range(T-1))
        if CAP_DDU:
            model.addConstrs(self.lowerR[r]-self.M*(1-self.e[t,r]) <= gu.quicksum(self.y[k] for k in range(t+1)) for r in range(R) for t in range(T))
            model.addConstrs(gu.quicksum(self.y[k] for k in range(t+1)) <= self.upperR[r]+self.M*(1-self.e[t,r])  for r in range(R) for t in range(T))

            model.addConstrs(gu.quicksum(self.e[t,r] for r in range(R))==1 for t in range(T) )
        
        
    def AddSupplyConstraints_Dual(self,model):
        #dual constraint 1
        model.addConstrs((-self.nu[t]+self.tau[t] <= self.c_VP[t] for t in range(T)), name='Dual_1')
        #dual constraints 3
        model.addConstrs((-self.tau[t]+self.outgoing[k].To.psi[t] <= self.c_VT[t]*self.outgoing[k].dist for k in range(demand_count) for t in range(T)), name='Dual_3' )

    def AddSupplyConstraints_Deterministic(self,model):
        #11b
        model.addConstrs(self.y[t]<=self.cap*self.x[t]  for t in range(T))
        #11c
        model.addConstrs(self.x[t]<=self.x[t+1] for t in range(T-1))
        
        model.addConstrs(self.h[t]<=  self.In_cap + gu.quicksum(self.y[k]/self.Efi[k] for k in range(t+1)) for t in range(T))
        #11e supply
        model.addConstrs(gu.quicksum(line.z[t] for line in self.outgoing) == self.h[t] for t in range(T))
        
    def AddSupplyConstraints_Simulation(self,model,Opt_y, s):
        
        for t in range(T):
            model.addConstr(self.h[s][t]<= self.In_cap + sum(Opt_y[self][k]/self.Efi[k] for k in range(t+1)), name="P_Cap"+str(self.node_id)+','+str(t)+','+str(s))
        #11e supply
        model.addConstrs(gu.quicksum(line.z[s][t] for line in self.outgoing) == self.h[s][t] for t in range(T))
        
    def AddDemandVariables_Budget(self, model, t):
        self.w[t] = model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="w["+str(self.node_id)+","+str(t)+"]")
        
 
    
    def AddPortConstraints_SP(self,model,t):
        #11e import
        model.addConstr(gu.quicksum(line.z[t] for line in self.outgoing) <= self.v[t]  )
        #dual constraints 2
        model.addConstr((self.ups[t] <= self.c_I[t] ),name='Dual_2')
        #dual constraints 4
        model.addConstrs((-self.ups[t]+self.outgoing[k].To.psi[t] <= self.outgoing[k].c_VT[t]*self.outgoing[k].dist for k in range(demand_count) ),name='Dual_4' )

        #complemantary slackness3
        model.addConstr((self.v[t]-gu.quicksum(line.z[t] for line in self.outgoing)<=self.M*self.CS_I1[t] ),name='CS_3')
        model.addConstr((self.ups[t]<=self.M*(1-self.CS_I1[t]) ), name='CS_3')
        
        #dual complementart slackness 2
        model.addConstr((self.c_I[t]-self.ups[t]<=self.M*self.D_I1[t] ),name='D_2')
        model.addConstr((self.v[t]<=self.M*(1-self.D_I1[t]) ),name='D_2')
        
        #dual complementart slackness 4
        model.addConstrs(self.outgoing[k].c_VT[t]*self.outgoing[k].dist+self.ups[t]-self.outgoing[k].To.psi[t]<=self.M*self.D_ID[t][k] for k in range(demand_count) )
        model.addConstrs(self.outgoing[k].z[t]<=self.M*(1-self.D_ID[t][k]) for k in range(demand_count) )

    def AddDemandConstraints_SP(self,model,t):
        model.addConstr((gu.quicksum(line.z[t] for line in self.incoming) >= self.xi[t] ), name='demand'+str(self.node_id)+','+str(t) )
        
        #complemantary slackness4
        model.addConstr((gu.quicksum(line.z[t] for line in self.incoming)-self.xi[t]<=self.M*self.CS_D1[t] ),name='CS_4' )
        model.addConstr((self.psi[t]<=self.M*(1-self.CS_D1[t]) ),name='CS_4')
    
    def AddPortConstraints_Dual(self,model):        
        #dual constraints 2
        model.addConstrs((self.ups[t] <= self.c_I[t] for t in range(T)),name='Dual_2')
        #dual constraints 4
        model.addConstrs((-self.ups[t]+self.outgoing[k].To.psi[t] <= self.c_VT[t]*self.outgoing[k].dist for k in range(demand_count) for t in range(T)),name='Dual_4' )

    def AddPortConstraints_Deterministic(self,model):
        model.addConstrs(gu.quicksum(line.z[t] for line in self.outgoing) <= self.v[t] for t in range(T)  )

    def AddPortConstraints_Simulation(self,model, s):
        model.addConstrs(gu.quicksum(line.z[s][t] for line in self.outgoing) <= self.v[s][t] for t in range(T))
        
    def AddDemandConstraints_MP(self,model):
        for line in self.incoming:
            if line.From in supply_nodes:
                if CAP_DDU:
                    model.addConstrs(self.zeta1_variables[line.From][t,r]<=self.beta_max*line.From.e[t,r] for r in range(R) for t in range(T))
                    model.addConstrs(self.zeta1_variables[line.From][t,r]<=self.beta1[t] for r in range(R) for t in range(T))
                    model.addConstrs(self.zeta1_variables[line.From][t,r]>=self.beta1[t] - self.beta_max*(1-line.From.e[t,r]) for r in range(R) for t in range(T))

                    model.addConstrs(self.zeta2_variables[line.From][t,r]<=self.beta_max*line.From.e[t,r] for r in range(R) for t in range(T))
                    model.addConstrs(self.zeta2_variables[line.From][t,r]<=self.beta2[t] for r in range(R) for t in range(T))
                    model.addConstrs(self.zeta2_variables[line.From][t,r]>=self.beta2[t] - self.beta_max*(1-line.From.e[t,r]) for r in range(R) for t in range(T))

                else:
                    model.addConstrs(self.phi1_variables[line.From][t]<=self.beta_max*line.From.x[t] for t in range(T))
                    model.addConstrs(self.phi1_variables[line.From][t]<=self.beta1[t] for t in range(T))
                    model.addConstrs(self.phi1_variables[line.From][t]>=self.beta1[t] - self.beta_max*(1-line.From.x[t])  for t in range(T))

                    model.addConstrs(self.phi2_variables[line.From][t]<=self.beta_max*line.From.x[t] for t in range(T))
                    model.addConstrs(self.phi2_variables[line.From][t]<=self.beta2[t] for t in range(T))
                    model.addConstrs(self.phi2_variables[line.From][t]>=self.beta2[t] - self.beta_max*(1-line.From.x[t])  for t in range(T))



def distribute_roles(supply_info, import_info, demand_info):  
   
    nodes = []; supply_nodes=[]; port_nodes=[]; demand_nodes=[]
    node_id=0
    for i in supply_info:
        nodes.append(Node(node_id, supply_info[i]['coords'][0], supply_info[i]['coords'][1], supply_info[i]['type'],0))
        supply_nodes.append(nodes[node_id])
        node_id+=1
    
    for i in import_info:
        nodes.append(Node(node_id, import_info[i]['coords'][0], import_info[i]['coords'][1], import_info[i]['type'],0))
        port_nodes.append(nodes[node_id])
        node_id+=1
        
    for i in demand_info:
        nodes.append(Node(node_id, demand_info[i]['coords'][0], demand_info[i]['coords'][1], demand_info[i]['type'],demand_info[i]['demand']))
        demand_nodes.append(nodes[node_id])
        node_id+=1
        
    return nodes,supply_nodes,port_nodes,demand_nodes

supply_count = 5
storage_count = 0
port_count =1
demand_count =13
total_locations = supply_count+storage_count+port_count+demand_count 
T=5
Hours=((2050-2030)/(T-1))*365*24
Int_Rate=0.95
PerPeriod=(2050-2030)/(T-1)
R= 3

supply_info = {
    0 : {'coords': [53.44059, 6.82363], 'type': 'production', 'descption': 'Eemshaven RWE Electrolyzer'},
    1 : {'coords': [53.319197, 6.944017], 'type': 'production', 'descption': 'Delfzijl Electrolyzer and By-product use'},
    2 : {'coords': [52.754560, 6.936359], 'type': 'production', 'descption': 'Emmen Electrolyzer'},
    3 : {'coords': [53.194191, 6.621600], 'type': 'production', 'descption': 'Groningen solar farm Electrolyzer'},
    4 : {'coords': [53.080842, 6.928074], 'type': 'production', 'descption': 'HyStock Electrolyzer'}
}

import_info = {
    0 : {'coords': [53.45162, 6.82991], 'type': 'port', 'descption': 'Eemshaven import terminal'},
    }

demand_info = {
    0  : {'coords': [53.316718, 6.952066], 'demand': [120000,260000], 'type': 'industry', 'descption': 'Industry cluster 1 in Delzijl (Chempark Delzijl) with RWE powerplant'},
    1  : {'coords': [53.44026, 6.83999],'demand': [50000,90000], 'type': 'industry', 'descption': 'Industry cluster 2 in Eemshaven '},
    2  : {'coords': [53.09345, 6.88291],'demand': [30000,60000], 'type': 'industry', 'descption': 'Industry cluster 3 in Veendam'},
    3  : {'coords': [52.774836, 6.906929],'demand': [150,400], 'type': 'industry', 'descption': 'Industry cluster 4 in Emmen'},
    4  : {'coords': [53.17392, 5.807429],'demand': [4000,12000], 'type': 'industry', 'descption': 'Industry cluster 5 in Leeuwarden'},
    5  : {'coords': [52.72928, 6.507629],'demand': [3000,9000], 'type': 'industry', 'descption': 'Industry cluster 6 in Hoogeveen'},
    6  : {'coords': [52.707468, 6.411857],'demand': [75,200], 'type': 'house', 'descption': 'Domestic Heating Hoogeveen'},
    7  : {'coords': [53.206287, 6.598420],'demand': [365,700], 'type': 'mobility', 'descption': 'Groningen Refueling Station/Groningen Depot'},
    8  : {'coords': [53.319264, 6.944039],'demand': [32.5,55], 'type': 'mobility', 'descption': 'Delfzijl Refueling Station'},
    9  : {'coords': [52.77276, 6.449729],'demand': [32.5,55], 'type': 'mobility', 'descption': 'Hoogeveen Refueling Station/Drenthe Depot'},
    10 : {'coords': [53.187391, 5.756540],'demand': [32.5,55], 'type': 'mobility', 'descption': 'Friesland Refueling Station Depot'},
    11 : {'coords': [53.12805, 6.58622],'demand': [16.43,27.38], 'type': 'airport', 'descption': 'Groningen Eelde airport'},
    12 : {'coords': [53.321321, 6.942968],'demand': [14.34,20], 'type': 'ship', 'descption': 'Delfzijl inland shipping refueling Station'},
}


nodes,supply_nodes,port_nodes,demand_nodes= distribute_roles(supply_info, import_info, demand_info)


opt_dro_ddu = {"y" : {supply_nodes[0]: {0: 0.0, 1: 0.0, 2: 231899.2165671924, 3: 163789.38820343284, 4: 0.0}, supply_nodes[1]: {0: 0.0, 1: 525645.5571070398, 2: 412694.850040042, 3: 187139.18048175084, 4: 217494.5979515799}, supply_nodes[2]: {0: 0.0, 1: 0.0, 2: 29068.734290291213, 3: 8694.9391824235, 4: 10676.733478166036}, supply_nodes[3]: {0: 0.0, 1: 24517.526502944937, 2: 15389.467010325194, 3: 11365.433256715136, 4: 14523.790238066667}, supply_nodes[4]: {0: 113013.698630137, 1: 85893.84496987921, 2: 34008.65849662308, 3: 39874.89893373153, 4: 45425.39261683514}},
               "x": {supply_nodes[0]: {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}, supply_nodes[1]: {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}, supply_nodes[2]: {0: 0.0, 1: 0.0, 2: 1.0, 3: 1.0, 4: 1.0}, supply_nodes[3]: {0: 0.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}, supply_nodes[4]: {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}}}
opt_dro= {"y": {supply_nodes[0]: {0: 0.0, 1: 0.0, 2: 262449.4728698448, 3: 5898.6343945728995, 4: 0.0}, supply_nodes[1]: {0: 0.0, 1: 525578.6727392139, 2: 107347.06626573663, 3: 124759.4536545005, 4: 184013.95561874233}, supply_nodes[2]: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}, supply_nodes[3]: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}, supply_nodes[4]: {0: 0.0, 1: 151873.9386668465, 2: 32068.160171291, 3: 75323.85906716734, 4: 68141.53047334531}},
          "x": {supply_nodes[0]: {0: 0.0, 1: 0.0, 2: 1.0, 3: 1.0, 4: 1.0}, supply_nodes[1]: {0: 0.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}, supply_nodes[2]: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}, supply_nodes[3]: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}, supply_nodes[4]: {0: 0.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}}}

class Line:
    def __init__(self, From, To):
        self.From=From
        self.To=To
        
        self.From_id=From.node_id
        self.To_id=To.node_id
  
        self.dist=Distance_Dict[From.node_id,To.node_id]

        From.outgoing.append(self)
        To.incoming.append(self)
        
        self.min=0
        self.max=float('inf')
        
        self.z={}
    
    def AddVariables_SP(self,model,t):
        #11h
        self.z[t]=model.addVar(lb=self.min,ub=self.max,vtype=gu.GRB.CONTINUOUS, name='z'+str(self.From_id)+','+str(self.To_id)+","+str(t) )
    
    def AddVariables_MP(self,model):
        self.MP_z={}

    def AddVariables_Dual(self,model):
        #11h
        self.z=model.addVars(T,lb=self.min,ub=self.max,vtype=gu.GRB.CONTINUOUS, name=['z'+str(self.From_id)+','+str(self.To_id)+","+str(t) for t in range(T) ])

    def AddVariables_Deterministic(self,model):
        self.z=model.addVars(T,lb=self.min,ub=self.max,vtype=gu.GRB.CONTINUOUS, name=['z'+str(self.From_id)+','+str(self.To_id)+","+str(t) for t in range(T) ])
    
    def AddVariables_Simulation(self,model, s):
        self.z[s]=model.addVars(T,lb=self.min,ub=self.max,vtype=gu.GRB.CONTINUOUS, name=['z'+str(self.From_id)+','+str(self.To_id)+","+str(t)+","+str(s) for t in range(T)])

Lines=[]
for i in nodes:
    for j in nodes:
        if i in supply_nodes or i in port_nodes:
            if j in demand_nodes:
                Lines.append(Line(i, j))

def AssignScenarios(demand_nodes, lambda_factor):
    for node in demand_nodes:
        factor=1000/365/24
        initial=node.Demand[0]*factor
        final=node.Demand[1]*factor
        alpha=(final/initial)**(1/(T-1))
        node.mu[0]= initial
        node.Support[0]=[initial*bound for bound in sup_range]
        node.epsilon[0]=initial*epsilon_factor
        for t in range(1,T):
            node.mu[t]=node.mu[t-1]*alpha
            node.Support[t]=[node.mu[t]*bound for bound in sup_range]
            node.epsilon[t]=node.mu[t]*epsilon_factor
        SUM=sum(math.exp(-line.dist / function)  for line in node.incoming  if line.From in supply_nodes and line.dist <=lambda_rate)
        if CAP_DDU:
            for line in node.incoming:
                if line.From in supply_nodes:
                    node.Lambda[line.From]={}
                    for xr in range(R):
                        if line.dist <=lambda_rate:
                            node.Lambda[line.From][xr]=(math.exp(-line.dist /function)/SUM)*lambda_factor*xr
                        else:
                            node.Lambda[line.From][xr]=0
        if not CAP_DDU:
            for line in node.incoming:
                if line.From in supply_nodes:
                    if line.dist <=lambda_rate:
                        node.Lambda[line.From]=(math.exp(-line.dist / function)/SUM)*lambda_factor
                    else:
                        node.Lambda[line.From]=0
        for node in supply_nodes:
            node.lowerR=[0,XR_R,2*XR_R]
            node.upperR=[XR_R,2*XR_R,node.M]
              
AssignScenarios(demand_nodes, lambda_factor)                  


def AssignCostParameters(supply_nodes,port_nodes,demand_nodes,Lines):
    for node in supply_nodes:
        node.c_FS={}  ; node.c_HS={}  ; node.c_VP={} ; node.Efi={}  
        node.c_FS[0]=capital ; node.c_HS[0]=1225
        node.c_VP[0]=3.8*Hours
        node.Efi[0]=33/0.8
        for t in range(1,T):
            node.c_FS[t]=node.c_FS[t-1]*(550/node.c_HS[0])**(1/(T-1))*(Int_Rate**PerPeriod)
            node.c_HS[t]=node.c_HS[t-1]*(550/node.c_HS[0])**(1/(T-1))*(Int_Rate**PerPeriod)
            node.c_VP[t]=node.c_VP[t-1]*((2.2*Hours)/node.c_VP[0])**(1/(T-1))*(Int_Rate**PerPeriod)
            node.Efi[t]=node.Efi[t-1]*((33/0.95)/node.Efi[0]) ** (1/(T-1))
    
    for node in port_nodes:   
        node.c_I={}
        node.c_I[0]=3.8*Hours
        for t in range(1,T):
            node.c_I[t]=node.c_I[t-1]*((2.2*C_IMPORT*Hours)/node.c_I[0]) ** (1/(T-1))*(Int_Rate**PerPeriod)
    
    for node in demand_nodes:
        node.R={}
        node.R[0]=profit*Hours
        for t in range(1,T):
            node.R[t]=node.R[t-1]*0.9*(Int_Rate**PerPeriod)
        
    for line in Lines:
        line.c_VT={}
        line.c_VT[0]=0.01*Hours
      
        for t in range(1,T):
            line.c_VT[t]=line.c_VT[t-1]*((0.008*Hours)/line.c_VT[0]) ** (1/(T-1))*(Int_Rate**PerPeriod)

            
AssignCostParameters(supply_nodes,port_nodes,demand_nodes,Lines)           

class MasterProblem:

    def __init__(self, opt_y = None, opt_x = None):
        self.fixed_sol = not (opt_y is None or opt_x is None)
        self.model= gu.Model()     
        self.model.setParam('OutputFlag', 0)
        self.model.setParam('Threads', multiprocessing.cpu_count())
        self.model.Params.MIPGap = 0.001

        for node in supply_nodes :
            node.AddSupplyVariables_MP(self.model, self.fixed_sol)

        for node in port_nodes :
            node.AddPortVariables_MP(self.model)  

        for node in demand_nodes :
            node.AddDemandVariables_MP(self.model, self.fixed_sol)

        for line in Lines :
            line.AddVariables_MP(self.model) 
        
        if not self.fixed_sol: 
                
            for node in supply_nodes :
                node.AddSupplyConstraints_MP(self.model)
            
            for node in demand_nodes :
                node.AddDemandConstraints_MP(self.model)
                

            self.opening=gu.quicksum(node.c_FS[0]*node.x[0] for node in supply_nodes )+ gu.quicksum(node.c_FS[t]*(node.x[t]-node.x[t-1]) for node in supply_nodes for t in range(1,T))
            self.capacity=gu.quicksum(node.c_HS[t]*node.y[t] for node in supply_nodes for t in range(T))
            if CAP_DDU:
                self.dual_beta1=gu.quicksum(node.mu[t]*node.beta1[t] +node.mu[t]*gu.quicksum(node.Lambda[supply][r]*node.zeta1_variables[supply][t,r] for supply in supply_nodes for r in range(R)) +node.epsilon[t]*node.beta1[t] for node in demand_nodes for t in range(T))
                self.dual_beta2=gu.quicksum(-node.mu[t]*node.beta2[t] -node.mu[t]*gu.quicksum(node.Lambda[supply][r]*node.zeta2_variables[supply][t,r] for supply in supply_nodes for r in range(R)) +node.epsilon[t]*node.beta2[t] for node in demand_nodes for t in range(T))
            else:
                self.dual_beta1=gu.quicksum(node.mu[t]*node.beta1[t] +node.mu[t]*gu.quicksum((node.Lambda[supply] if DDU else 0)*node.phi1_variables[supply][t] for supply in supply_nodes) +node.epsilon[t]*node.beta1[t] for node in demand_nodes for t in range(T))
                self.dual_beta2=gu.quicksum(-node.mu[t]*node.beta2[t] -node.mu[t]*gu.quicksum((node.Lambda[supply] if DDU else 0)*node.phi2_variables[supply][t] for supply in supply_nodes) +node.epsilon[t]*node.beta2[t] for node in demand_nodes for t in range(T))
        else:
            self.Opt_y = opt_y
            self.Opt_x = opt_x

            self.opening=gu.quicksum(node.c_FS[0]*opt_x[node][0] for node in supply_nodes )+ gu.quicksum(node.c_FS[t]*(opt_x[node][t]-opt_x[node][t-1]) for node in supply_nodes for t in range(1,T))
            self.capacity=gu.quicksum(node.c_HS[t]*opt_y[node][t] for node in supply_nodes for t in range(T))
            self.dual_beta1=gu.quicksum(node.mu[t]*node.beta1[t] +node.mu[t]*node.beta1[t]*gu.quicksum(node.Lambda[supply]*opt_x[supply][t] for supply in supply_nodes) +node.epsilon[t]*node.beta1[t] for node in demand_nodes for t in range(T))
            self.dual_beta2=gu.quicksum(-node.mu[t]*node.beta2[t] -node.mu[t]*node.beta2[t]*gu.quicksum(node.Lambda[supply]*opt_x[supply][t] for supply in supply_nodes) +node.epsilon[t]*node.beta2[t] for node in demand_nodes for t in range(T))
        
        self.alpha=self.model.addVar(lb=-gu.GRB.INFINITY, ub=gu.GRB.INFINITY, vtype=gu.GRB.CONTINUOUS,name="alpha")

        self.total= self.opening+self.capacity + self.alpha +self.dual_beta1+ self.dual_beta2
    
        self.model.setObjective( self.total ,gu.GRB.MINIMIZE)

                        
    def optimize(self,time_limit):
        self.model.setParam('TimeLimit', time_limit)
        self.start = timeit.default_timer()
        self.model.optimize()
        self.stop = timeit.default_timer()
        self.GurobiTime=self.stop-self.start
        self.model.write(f"MasterModel_{int(self.fixed_sol)}.lp")
        
        if self.model.status == gu.GRB.TIME_LIMIT:
            return self.Opt_y,   self.Opt_beta1,  self.Opt_beta2 
            
        self.obj=self.model.objval
        self.bound=self.model.objBound
        
        if not self.fixed_sol:
            self.Opt_x={}
            self.Opt_y={}
            
            for node in supply_nodes:
                self.Opt_y[node]={}
                self.Opt_x[node]={}
                for t in range(T):
                    self.Opt_y[node][t]=node.y[t].X
                    self.Opt_x[node][t]=node.x[t].X
                
        self.Opt_beta1={}
        self.Opt_beta2={}
        for node in demand_nodes:
            self.Opt_beta1[node]={}
            self.Opt_beta2[node]={}
            for t in range(T):
                self.Opt_beta1[node][t]=node.beta1[t].X
                self.Opt_beta2[node][t]=node.beta2[t].X
                
        self.SC=self.model.SolCount 
        
        self.opt_alpha=self.alpha.X
        return  self.Opt_y,   self.Opt_beta1,  self.Opt_beta2 
            
        
    def Start(self):
        for node in supply_nodes:
            for t in range(T):
                node.y[t].start=self.Opt_y[node][t]
        
        for node in demand_nodes:
            for t in range(T):
                node.beta1[t].start=self.Opt_beta1[node][t]
                node.beta2[t].start=self.Opt_beta2[node][t]
        
    def AddColumnandConstraints(self,It,XI):
        for node in supply_nodes :
            node.MP_h[It]=self.model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["h"+str(node.node_id)+","+str(t) for t in range(T)])    
                        
        for node in port_nodes :
            node.MP_v[It]= self.model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["v"+str(node.node_id)+","+str(t)  for t in range(T)])
            
        for line in Lines :
            line.MP_z[It]=self.model.addVars(T,lb=line.min,ub=line.max,vtype=gu.GRB.CONTINUOUS, name=['z'+str(line.From_id)+','+str(line.To_id)+","+str(t) for t in range(T) ])

        for node in supply_nodes:
            for t in range(T):
                if self.fixed_sol:
                    self.model.addConstr((gu.quicksum(self.Opt_y[node][k]/node.Efi[k] for k in range(t+1)) +node.In_cap >=node.MP_h[It][t]) , name="P_Cap"+str(It)+','+str(node.node_id)+','+str(t))
                else:
                    self.model.addConstr((gu.quicksum(node.y[k]/node.Efi[k] for k in range(t+1)) +node.In_cap >=node.MP_h[It][t]) , name="P_Cap"+str(It)+','+str(node.node_id)+','+str(t))
                self.model.addConstr((gu.quicksum(line.MP_z[It][t] for line in node.outgoing) <= node.MP_h[It][t]), name="Production"+str(It)+','+str(node.node_id)+','+str(t) )
        
        for node in port_nodes :
            for t in range(T):
                self.model.addConstr((gu.quicksum(line.MP_z[It][t] for line in node.outgoing) <= node.MP_v[It][t] ), name="Import"+str(It)+','+str(node.node_id)+','+str(t) )

        for node in demand_nodes :
            for t in range(T):
                self.model.addConstr((gu.quicksum(line.MP_z[It][t] for line in node.incoming) >= XI[node][t] ), name="Demad"+str(It)+','+str(node.node_id)+','+str(t) )


        self.model.addConstr(self.alpha+gu.quicksum((node.beta1[t]-node.beta2[t])*XI[node][t] for t in range(T) for node in demand_nodes) >= gu.quicksum(node.c_VP[t]*node.MP_h[It][t] for node in supply_nodes for t in range(T)) + gu.quicksum(node.c_I[t]*node.MP_v[It][t] for node in port_nodes for t in range(T) ) + gu.quicksum(-node.R[t]* XI[node][t] for node in demand_nodes for t in range(T)) + gu.quicksum(line.c_VT[t]*line.dist*line.MP_z[It][t] for line in Lines for t in range(T)), name='cut'+str(It) )
        
        self.model.update()

    def BendersCut(self,It,XI,SP_nu,SP_psi):
        self.model.addConstr(self.alpha+gu.quicksum((node.beta1[t]-node.beta2[t])*XI[node][t] for t in range(T) for node in demand_nodes) >= -gu.quicksum( gu.quicksum( node.y[k] for k in range(t+1))* SP_nu[node][t] for node in supply_nodes for t in range(T)) + gu.quicksum( XI[node][t]*(SP_psi[node][t]-node.R[t]) for node in demand_nodes for t in range(T)) ) 
    
    def UpdateObj(self):
        self.opening=gu.quicksum(node.c_FS[0]*node.x[0] for node in supply_nodes )+ gu.quicksum(node.c_FS[t]*(node.x[t]-node.x[t-1]) for node in supply_nodes for t in range(1,T))
        self.capacity=gu.quicksum(node.c_HS[t]*node.y[t] for node in supply_nodes for t in range(T))
        self.dual_beta1=gu.quicksum(node.mu[t]*node.beta1[t] +node.mu[t]*gu.quicksum(node.Lambda[supply][r]*node.zeta1_variables[supply][t,r] for supply in supply_nodes for r in range(R)) +node.epsilon[t]*node.beta1[t] for node in demand_nodes for t in range(T))
        self.dual_beta2=gu.quicksum(-node.mu[t]*node.beta2[t] -node.mu[t]*gu.quicksum(node.Lambda[supply][r]*node.zeta2_variables[supply][t,r] for supply in supply_nodes for r in range(R)) +node.epsilon[t]*node.beta2[t] for node in demand_nodes for t in range(T))
        self.total= self.opening+self.capacity + self.alpha +self.dual_beta1+ self.dual_beta2    
        self.model.setObjective( self.total ,gu.GRB.MINIMIZE)

class BendersSubProblem:

    def __init__(self,Lindex):
        self.t=Lindex
        self.model= gu.Model()  
        self.model.setParam('OutputFlag', 0)
        self.model.Params.MIPGap = 0.0001 

        for node in supply_nodes :
            node.AddSupplyVariables_SP(self.model,self.t)
                        
        for node in port_nodes :
            node.AddPortVariables_SP(self.model,self.t)
        
        for node in demand_nodes :
            node.AddDemandVariables_SP(self.model,self.t)
            
        for line in Lines :
            line.AddVariables_SP(self.model,self.t)
        
        for node in supply_nodes :
            node.AddSupplyConstraints_SP(self.model,self.t)
        
        for node in port_nodes:
            node.AddPortConstraints_SP(self.model,self.t)
        
        for node in demand_nodes:
            node.AddDemandConstraints_SP(self.model,self.t)
        
        self.productioncost= gu.quicksum(node.c_VP[self.t]*node.h[self.t] for node in supply_nodes )
        self.importcost= gu.quicksum(node.c_I[self.t]*node.v[self.t] for node in port_nodes )
        self.transportationcost=gu.quicksum(line.c_VT[self.t]*line.dist*line.z[self.t] for line in Lines )
        self.betasum=gu.quicksum(-(Opt_beta1[node][self.t]-Opt_beta2[node][self.t])*node.xi[self.t] for node in demand_nodes )
        self.totalreward= gu.quicksum(-node.R[self.t]* node.xi[self.t] for node in demand_nodes ) 
        
        self.total=self.productioncost+self.importcost+self.transportationcost+self.betasum+self.totalreward
        self.model.setObjective(self.total,gu.GRB.MAXIMIZE)
                                         
        self.model.update()

    def optimize(self,time_limit):
        self.model.setParam('TimeLimit', time_limit)
        self.start = timeit.default_timer()
        self.model.optimize()
        self.stop = timeit.default_timer()
        self.GurobiTime=self.stop-self.start
        self.model.write("1SPmodel.lp")
        
        if self.model.status == 3 or self.model.status == 4 :
            self.model.computeIIS()
            self.model.write("1SPmodel.ilp")
            
        
        self.productioncost_opt=self.productioncost.getValue()
        self.importcost_opt=self.importcost.getValue()
        self.transportationcost_opt=self.transportationcost.getValue()
        self.betasum_opt=self.betasum.getValue()
        self.totalreward_opt=self.totalreward.getValue()
        self.total_opt=self.total.getValue()
                
        self.SC=self.model.SolCount 
        
        self.obj={}
        
        
        self.XI={}
        self.opt_nu={}
        self.opt_psi={}
        for k in range(self.SC):
            self.model.Params.SolutionNumber = k
            self.obj[k]=self.model.PoolObjVal
            self.XI[k]={}
            for node in demand_nodes:
                self.XI[k][node]=node.xi[self.t].xn
            self.opt_nu[k]={}  
            for node in supply_nodes:
                self.opt_nu[k][node]=node.nu[self.t].xn
            self.opt_psi[k]={}
            for node in demand_nodes:
                self.opt_psi[k][node]=node.psi[self.t].xn

  
    def ChangeRHS(self,Opt_y):
        for node in supply_nodes:
            self.model.getConstrByName("P_Cap"+str(node.node_id)+','+str(self.t)).rhs= -(sum(Opt_y[node][k]/node.Efi[k] for k in range(self.t+1)) +node.In_cap)
            self.model.getConstrByName('CS_S1C'+str(node.node_id)+','+str(self.t)).rhs= -(sum(Opt_y[node][k]/node.Efi[k] for k in range(self.t+1)) +node.In_cap)
    
    def ChangeObjective(self,Opt_beta1, Opt_beta2):
        for node in demand_nodes:
            node.xi[self.t].Obj=-(Opt_beta1[node][self.t]-Opt_beta2[node][self.t])-node.R[self.t]

def SPsolve(t: int, time_limit: int) -> int:
    SP[t].ChangeRHS(Opt_y)
    SP[t].ChangeObjective(Opt_beta1, Opt_beta2)
    SP[t].optimize(time_limit) 
    return 0

Mstart = timeit.default_timer()
MP=MasterProblem()

Warm_XI={}
for node in demand_nodes:
    Warm_XI[node]={}
    for t in range(T):
        Warm_XI[node][t]=node.Support[t][-1]        

Warm_XI2={}
for node in demand_nodes:
    Warm_XI2[node]={}
    for t in range(T):
        Warm_XI2[node][t]=node.Support[t][0]

MP.AddColumnandConstraints(-1,Warm_XI)
MP.AddColumnandConstraints(-2,Warm_XI2)

WARM={}
for wit in range(len(demand_nodes)*T):
    WARM[wit]={}
    for node in demand_nodes:
        WARM[wit][node]={}
        for t in range(T):
            WARM[wit][node][t]=node.Support[t][0]

wit=0
for i in range(len(demand_nodes)):
    for t in range(T):
        WARM[wit][demand_nodes[i]][t]=demand_nodes[i].Support[t][-1]
        wit+=1
  

Opt_y, Opt_beta1, Opt_beta2 = MP.optimize(TIME_LIMIT)  
Mstop = timeit.default_timer()

SP=[BendersSubProblem(t)for t in range(T)]

UB = float("inf")
LB = -float("inf")
AIt=0
It=0
perc5 = float("inf") ; perc2= float("inf") ; perc1=float("inf")
time5=float("inf") ;time2=float("inf") ; time1=float("inf") 
In_tol=0.50
End_tol=0.001

start = timeit.default_timer()
while (UB-LB)>=abs(0.001*UB):
    MP.model.Params.MIPGap = max(min(In_tol,0.50*(UB-LB)/abs(UB)),End_tol) 
    print('gap',max(min(In_tol,0.50*(UB-LB)/abs(UB)),End_tol))
    if It!=0:
        Opt_y, Opt_beta1, Opt_beta2=MP.optimize(TIME_LIMIT - (timeit.default_timer()-start)) 
        MP.Start()
    LB=MP.bound
    SPstart = timeit.default_timer()    
    Opt_XI={}
    time_limit = TIME_LIMIT - (timeit.default_timer()-start)
    if time_limit < 0:
        break
    if BUDGET == 100:
        for t in range(T):
            SPsolve(t, time_limit)
        SPstop = timeit.default_timer()  

        SP_sum=sum(SP[t].obj[0] for t in range(T))
        if SP_sum+(MP.obj-MP.opt_alpha)<=UB:
                UB=SP_sum+(MP.obj-MP.opt_alpha)
                optimalY=Opt_y    
                optimalX=MP.Opt_x
                optimalBeta1=Opt_beta1
                optimalBeta2=Opt_beta2
                optimalAlpha = MP.opt_alpha
        for node in demand_nodes:
            Opt_XI[node]={}
            for t in range(T):
                Opt_XI[node][t]=SP[t].XI[0][node]
    else:
        SP.ChangeRHS(Opt_y)
        SP.ChangeObjective(Opt_beta1, Opt_beta2)
        SP.optimize(time_limit) 
        SPstop = timeit.default_timer()  
        SP_sum=SP.obj[0]
        if SP_sum+(MP.obj-MP.opt_alpha)<=UB:
                UB=SP_sum+(MP.obj-MP.opt_alpha)
                optimalY=Opt_y
                optimalX=MP.Opt_x
                optimalBeta1=MP.Opt_beta1
                optimalBeta2=MP.Opt_beta2
        for node in demand_nodes:
            Opt_XI[node]={}
            for t in range(T):
                Opt_XI[node][t]=SP.XI[0][t][node]
    MP.AddColumnandConstraints(It,Opt_XI)    
    It+=1
    if 100*(UB-LB)/abs(UB)<= 5:
        perc5=min(It,perc5)
        time5=min(timeit.default_timer()-start,time5)
    if 100*(UB-LB)/abs(UB)<= 2:
        perc2=min(It,perc2)
        time2=min(timeit.default_timer()-start,time2)
    if 100*(UB-LB)/abs(UB)<= 1:
        perc1=min(It,perc1)
        time1=min(timeit.default_timer()-start,time1)
    print(It,LB,UB,100*(UB-LB)/abs(UB))

stop = timeit.default_timer()

print('time',stop-start)

dro_opt_y = {node.node_id: [v for v in val.values()] for node, val in optimalY.items()}
dro_opt_x = {node.node_id: [v for v in val.values()] for node, val in optimalX.items()}


out=[sum(optimalY[node][t]/1000 for node in supply_nodes) for t in range(T)]
List=[XR_R,LB,UB,It,stop-start,perc5,time5,perc2,time2,perc1,time1,dro_opt_y,dro_opt_x,out]
with open('simulation_csp_ddu.csv', 'a') as f_object: 
    writer_object = writer(f_object)
    writer_object.writerow(List)
    f_object.close()
