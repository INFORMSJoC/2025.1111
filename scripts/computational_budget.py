import numpy as np
import gurobipy as gu
import math
import timeit
import math
import argparse
from csv import writer
from enum import Enum

TIME_LIMIT = 18000

class Mode(Enum):
    CCG = 0
    CCGPlus = 1
    BendersDecomposition = 2
    ReformulatedSubProblem = 3
    BilinearSubProblem = 4
    ReformulatedSubProblemPlus = 5
    BilinearSubProblemPlus = 6

parser = argparse.ArgumentParser(description='Job number')
parser.add_argument('--job', '-d', metavar='the-no-of-data', dest='job_num', type=int, help='Choose job number', default=0)
args = parser.parse_args()

JOB=args.job_num

JOBlist=[]
for s in range(4,11, 2):
    for t in range(3,6):
        for b in [20,40,60,80]:
            for j in [2,3,4]:
                JOBlist.append([s,t,b,j])

supply_count=JOBlist[JOB][0]
demand_count=supply_count*2
T=JOBlist[JOB][1]
mode=Mode.CCG
budget=JOBlist[JOB][2]/100
rand_seed=JOBlist[JOB][3]

class Node:
    def __init__(self, node_id, x, y, role):
        self.node_id = node_id
        self.Coord_x = x
        self.Coord_y = y
        self.role = role
        self.outgoing =[]
        self.incoming=[]
        self.cap=10000     
        self.hd={}
        self.hprop={}
        #initial capacity of the supply node at the beginning
        self.In_cap=0
        self.w_min=0
        self.w_max=float('inf')
        self.w_0=0
        
        self.P_vec={}
        
        self.beta_max=10000
        self.psi_max=10000
        self.Support={}
        self.mu={}
        self.Lambda={}
        self.epsilon={}
        self.M=10000
        
        self.h={} ; self.nu={} ; self.tau={} ; self.CS_S1={} ;  self.CS_S2={} ;  self.D_S1={} ;  self.D_SD={}
        
        self.v={} ; self.ups={} ; self.CS_I1={} ; self.D_I1= {} ; self.D_ID={}

        self.xi={} ; self.psi={} ; self.CS_D1={} ; self.b0 = {} ; self.gamma = {} ; self.chi={} ; self.w={}

 
         
    
    def AddSupplyVariables_SP(self,model,t):
        self.h[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="h"+str(self.node_id)+","+str(t) )    
        self.nu[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="nu"+str(self.node_id)+","+str(t))    
        self.tau[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="tau"+str(self.node_id)+","+str(t) )    
    
        self.CS_S1[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="CS_S1"+str(self.node_id)+","+str(t) )    
        self.CS_S2[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="CS_S2"+str(self.node_id)+","+str(t) )    

        self.D_S1[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="D_S1"+str(self.node_id)+","+str(t) )            
        self.D_SD[t]=model.addVars(demand_count,lb=0,vtype=gu.GRB.BINARY,name=["D_SD"+str(self.node_id)+","+str(t)+","+str(k) for k in range(demand_count) ])    


    def AddSupplyVariables_MP(self,model):
        self.x=model.addVars(T,vtype=gu.GRB.BINARY,name=["x"+str(self.node_id)+","+str(t) for t in range(T)])
        self.y=model.addVars(T,lb=0,ub=self.cap,vtype=gu.GRB.CONTINUOUS,name=["y"+str(self.node_id)+","+str(t) for t in range(T)])

        self.MP_h={}
    
    def AddSupplyVariables_Dual(self,model, t2=-1):
        if t2 == -1:
            self.nu=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["nu"+str(self.node_id)+","+str(t) for t in range(T)])    
            self.tau=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["tau"+str(self.node_id)+","+str(t) for t in range(T)])    
        else:
            self.nu[t2]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="nu"+str(self.node_id)+","+str(t2))    
            self.tau[t2]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="tau"+str(self.node_id)+","+str(t2) )  

                  
    def AddPortVariables_SP(self,model,t):
        self.v[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="v"+str(self.node_id)+","+str(t)  )
        self.ups[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="ups"+str(self.node_id)+","+str(t)  )
        
        self.CS_I1[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="CS_I1"+str(self.node_id)+","+str(t) )    
        self.D_I1[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="D_I1"+str(self.node_id)+","+str(t) )    
        self.D_ID[t]=model.addVars(demand_count,lb=0,vtype=gu.GRB.BINARY,name=["D_I1"+str(self.node_id)+","+str(t)+","+str(k) for k in range(demand_count) ]) 

    def AddPortVariables_MP(self,model):
        self.MP_v={}
        
    def AddPortVariables_Dual(self,model, t2 = -1):
        if t2 == -1:
            self.ups=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["ups"+str(self.node_id)+","+str(t)  for t in range(T)])
        else:
            self.ups[t2]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="ups"+str(self.node_id)+","+str(t2))

                
    def AddDemandVariables_MP(self,model):
        self.beta1=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["beta1"+","+str(self.node_id)+","+str(t) for t in range(T) ])
        self.beta2=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["beta2"+","+str(self.node_id)+","+str(t) for t in range(T) ])
        self.phi1_variables={}
        self.phi2_variables={}
        for line in self.incoming:
            if line.From in supply_nodes:
                self.phi1_variables[line.From] = model.addVars(T,lb=0, vtype=gu.GRB.CONTINUOUS, name= ["phi1"+","+str(self.node_id)+","+str(t) for t in range(T) ])
                self.phi2_variables[line.From] = model.addVars(T,lb=0, vtype=gu.GRB.CONTINUOUS,  name= ["phi2"+","+str(self.node_id)+","+str(t) for t in range(T) ])

    def AddDemandVariables_SP(self,model,t):
        self.xi[t]=model.addVar(lb=self.Support[t][0],ub=self.Support[t][-1], vtype=gu.GRB.CONTINUOUS, name="xi"+","+str(self.node_id)+","+str(t))
        self.psi[t]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="psi"+str(self.node_id)+","+str(t) )

        self.CS_D1[t]=model.addVar(lb=0,vtype=gu.GRB.BINARY,name="CS_D1"+str(self.node_id)+","+str(t))    
        
    def AddDemandVariables_Dual(self,model, t2 = -1):
        if t2==-1:
            for t in range(T):
                self.xi[t]=model.addVar(lb=self.Support[t][0],ub=self.Support[t][-1], vtype=gu.GRB.CONTINUOUS,name="xi"+","+str(self.node_id)+","+str(t))

            self.psi=model.addVars(T,lb=0,vtype=gu.GRB.CONTINUOUS,name=["psi"+str(self.node_id)+","+str(t)  for t in range(T)])
        else:
            self.xi[t2]=model.addVar(lb=self.Support[t2][0],ub=self.Support[t2][-1], vtype=gu.GRB.CONTINUOUS,name="xi"+","+str(self.node_id)+","+str(t2))
            self.psi[t2]=model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="psi"+str(self.node_id)+","+str(t2))
        
    def AddDemandVariables_Alt(self, model, t):
        self.b0[t] = model.addVar(vtype=gu.GRB.BINARY,name="b0["+str(self.node_id)+","+str(t)+"]")
        self.gamma[t] = model.addVar(lb=0,ub=self.psi_max,vtype=gu.GRB.CONTINUOUS,name="gamma["+str(self.node_id)+","+str(t)+"]")
        self.chi[t] = self.R[t]

    def AddDemandVariables_Budget(self, model, t):
        self.w[t] = model.addVar(lb=0,vtype=gu.GRB.CONTINUOUS,name="w["+str(self.node_id)+","+str(t)+"]")
        
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
        
        
    def AddSupplyConstraints_Dual(self,model, t2 = -1):
        #dual constraint 1
        if t2 == -1:
            model.addConstrs((-self.nu[t]+self.tau[t] <= self.c_VP[t] for t in range(T)), name='Dual_1')
            #dual constraints 3
            model.addConstrs((-self.tau[t]+self.outgoing[k].To.psi[t] <= self.outgoing[k].c_VT[t]*self.outgoing[k].dist for k in range(demand_count) for t in range(T)), name='Dual_3' )
        else:
            model.addConstr(-self.nu[t2]+self.tau[t2] <= self.c_VP[t2], name='Dual_1')
            #dual constraints 3
            model.addConstrs((-self.tau[t2]+self.outgoing[k].To.psi[t2] <= self.outgoing[k].c_VT[t2]*self.outgoing[k].dist for k in range(demand_count) ), name='Dual_3' )

        
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
    
    def AddPortConstraints_Dual(self,model, t2 = -1):        
        #dual constraints 2
        if t2==-1:
            model.addConstrs((self.ups[t] <= self.c_I[t] for t in range(T)),name='Dual_2')
            #dual constraints 4
            model.addConstrs((-self.ups[t]+self.outgoing[k].To.psi[t] <= self.outgoing[k].c_VT[t]*self.outgoing[k].dist for k in range(demand_count) for t in range(T)),name='Dual_4' )
        else:
            model.addConstr((self.ups[t2] <= self.c_I[t2]),name='Dual_2')
            #dual constraints 4
            model.addConstrs((-self.ups[t2]+self.outgoing[k].To.psi[t2] <= self.outgoing[k].c_VT[t2]*self.outgoing[k].dist for k in range(demand_count)),name='Dual_4' )
    
    def AddDemandConstraints_MP(self,model):
        for line in self.incoming:
            if line.From in supply_nodes:
                model.addConstrs(self.phi1_variables[line.From][t]<=self.beta_max*line.From.x[t] for t in range(T))
                model.addConstrs(self.phi1_variables[line.From][t]<=self.beta1[t] for t in range(T))
                model.addConstrs(self.phi1_variables[line.From][t]>=self.beta1[t] - self.beta_max*(1-line.From.x[t])  for t in range(T))

                model.addConstrs(self.phi2_variables[line.From][t]<=self.beta_max*line.From.x[t] for t in range(T))
                model.addConstrs(self.phi2_variables[line.From][t]<=self.beta2[t] for t in range(T))
                model.addConstrs(self.phi2_variables[line.From][t]>=self.beta2[t] - self.beta_max*(1-line.From.x[t])  for t in range(T))

    def AddDemandConstraints_Alt(self,model, t):
        model.addConstr(self.psi[t] - self.M*self.b0[t] >= self.chi[t] - self.M, name = f'Psi[{self.node_id},{t}]')
        model.addConstr(self.gamma[t] <= self.psi_max*self.b0[t])
        model.addConstr(self.gamma[t] <= self.psi[t])
        model.addConstr(self.gamma[t] >= self.psi[t] - self.psi_max*(1 - self.b0[t]))
        
    def AddDemandConstraints_Budget(self, model, t):
        model.addConstr(self.w[t] >= -(self.xi[t] - (self.Support[t][0] + self.Support[t][-1]) /2))
        model.addConstr(self.w[t] >= self.xi[t] - (self.Support[t][0] + self.Support[t][-1]) /2)

     
def distribute_roles(total_locations, supply_locations, storage_locations, port_locations, demand_locations):  
    # Create a list of available roles based on the given counts
    target_sum = supply_locations+storage_locations+port_locations+demand_locations
    elements = [np.random.randint(1, 5) for _ in range(total_locations)]

    while sum(elements) != target_sum:
        index =  np.random.randint(0, total_locations)
        if sum(elements) > target_sum and elements[index] > 1:
            elements[index] -= 1
        elif sum(elements) < target_sum:
            elements[index] += 1
    elements=sorted(elements, reverse=True)
    
    nodes = []
    
    available_supply=["Supply"]*supply_locations
    available_storage=["Storage"]*storage_locations
    available_port=["Port"]*port_locations
    available_demand=["Demand"]*demand_locations
    
    for node_id in range(0, total_locations):               
        x = np.random.uniform(0, 100)   # Adjust the coordinate range as needed
        y = np.random.uniform(0, 100) 
        
        available_roles=[]
        if available_supply:
            available_roles.append("Supply")
        if available_storage:
            available_roles.append("Storage")
        if available_port: 
            available_roles.append("Port")
        if available_demand:
            available_roles.append("Demand")
            
        if available_roles:
            role = np.random.choice(available_roles, size=elements[node_id], replace=False) # Assign one available role to each node
            #to create a node for each role: for r in role:
            nodes.append(Node(node_id, x, y, role))
        
        if "Supply" in role:
            available_supply.remove("Supply")        
        if "Storage" in role:
            available_storage.remove("Storage")       
        if "Port" in role:
            available_port.remove("Port")                         
        if "Demand" in role:
            available_demand.remove("Demand")

    
    return nodes,available_roles


np.random.seed(rand_seed)

port_count=max(1,round(supply_count/10))
storage_count = 0
total_locations = supply_count + storage_count + port_count + demand_count 

nodes,av = distribute_roles(total_locations, supply_count, storage_count, port_count, demand_count)


supply_nodes=[]; storage_nodes=[]; port_nodes=[]; demand_nodes=[]
for node in nodes:
    if "Supply" in node.role:
        supply_nodes.append(node)
    if "Storage" in node.role:
        storage_nodes.append(node)
    if "Port" in node.role:
        port_nodes.append(node)                      
    if "Demand" in node.role:
        demand_nodes.append(node)

class Line:
    def __init__(self, From, To):
        self.From=From
        self.To=To
        
        self.From_id=From.node_id
        self.To_id=To.node_id
                
        self.dist= math.sqrt((From.Coord_x-To.Coord_x)**2 + (From.Coord_y-To.Coord_y)**2)

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

Lines=[]
for i in nodes:
    for j in nodes:
        if i in supply_nodes or i in port_nodes or i in storage_nodes:
            if j in demand_nodes:
                Lines.append(Line(i, j))
        if i in supply_nodes or i in port_nodes:
            if j in storage_nodes:
                Lines.append(Line(i, j))

def AssignScenarios(demand_nodes):
    for node in demand_nodes:        
        alpha=1.5
        initial=np.random.normal(50,5,1)[0]
        
        node.mu[0]= initial
        node.Support[0]=[initial*0.75,initial*1.25]
        node.epsilon[0]=initial*0.5
        for t in range(1,T):
            node.mu[t]=round(node.mu[t-1]*alpha,1)
            node.Support[t]=[node.mu[t]*0.75,node.mu[t]*1.25]
            node.epsilon[t]=node.mu[t]*0.5
        SUM=sum(math.exp(-line.dist / 25)  for line in node.incoming  if line.From in supply_nodes)
        for line in node.incoming:
            if line.From in supply_nodes:
                node.Lambda[line.From]=(math.exp(-line.dist / 25)/SUM)/2
              
AssignScenarios(demand_nodes)                  


def AssignCostParameters(supply_nodes,port_nodes,demand_nodes,Lines):
    for node in supply_nodes:
        node.c_FS={}  ; node.c_HS={}  ; node.c_VP={}        
        node.c_FS[0]=np.random.normal(4000,400,1)[0] ; node.c_HS[0]=np.random.normal(20,2,1)[0] ; node.c_VP[0]=10
        for t in range(1,T):
            node.c_FS[t]=node.c_FS[t-1]*0.5
            node.c_HS[t]=node.c_HS[t-1]*0.5
            node.c_VP[t]=node.c_VP[t-1]*0.9
    
    for node in port_nodes:   
        node.c_I={}
        node.c_I[0]=150
        for t in range(1,T):
            node.c_I[t]=node.c_I[t-1]*0.9
    
    for node in demand_nodes:
        node.R={}
        node.R[0]=100
        for t in range(1,T):
            node.R[t]=node.R[t-1]*0.9
        
    for line in Lines:
        line.c_VT={}
        line.c_VT[0]=0.5
        for t in range(1,T):
            line.c_VT[t]=line.c_VT[t-1]*0.9

AssignCostParameters(supply_nodes,port_nodes,demand_nodes,Lines)

class MasterProblem:

    def __init__(self):
        self.model= gu.Model()     
        self.model.setParam('OutputFlag', 0)

        self.model.Params.MIPGap = 0.001
        self.model.setParam('TimeLimit', 10000)
        
        for node in supply_nodes :
            node.AddSupplyVariables_MP(self.model)
                        
        for node in port_nodes :
            node.AddPortVariables_MP(self.model)
            
        for node in demand_nodes :
            node.AddDemandVariables_MP(self.model)
            
        for line in Lines :
            line.AddVariables_MP(self.model)   
            
        for node in supply_nodes :
            node.AddSupplyConstraints_MP(self.model)
        
        for node in demand_nodes :
            node.AddDemandConstraints_MP(self.model)
            
        self.alpha=self.model.addVar(lb=-gu.GRB.INFINITY, ub=gu.GRB.INFINITY, vtype=gu.GRB.CONTINUOUS,name="alpha")
       
        self.opening=gu.quicksum(node.c_FS[0]*node.x[0] for node in supply_nodes )+ gu.quicksum(node.c_FS[t]*(node.x[t]-node.x[t-1]) for node in supply_nodes for t in range(1,T))
        self.capacity=gu.quicksum(node.c_HS[t]*node.y[t] for node in supply_nodes for t in range(T))
        self.dual_beta1=gu.quicksum(node.mu[t]*node.beta1[t] +node.mu[t]*gu.quicksum(node.Lambda[supply]*node.phi1_variables[supply][t] for supply in supply_nodes) +node.epsilon[t]*node.beta1[t] for node in demand_nodes for t in range(T))
        self.dual_beta2=gu.quicksum(-node.mu[t]*node.beta2[t] -node.mu[t]*gu.quicksum(node.Lambda[supply]*node.phi2_variables[supply][t] for supply in supply_nodes) +node.epsilon[t]*node.beta2[t] for node in demand_nodes for t in range(T))

        self.total= self.opening+self.capacity + self.alpha +self.dual_beta1+ self.dual_beta2
    
        self.model.setObjective( self.total ,gu.GRB.MINIMIZE)
                        
    def optimize(self,time_limit):
        self.model.setParam('TimeLimit', time_limit)

        self.start = timeit.default_timer()
        self.model.optimize()
        self.stop = timeit.default_timer()
        self.GurobiTime=self.stop-self.start
        self.model.write("MasterModel.lp")

        if self.model.status == gu.GRB.TIME_LIMIT:
            return self.Opt_y,   self.Opt_beta1,  self.Opt_beta2 

        if self.model.status == 3 or self.model.status == 4 :
            self.model.computeIIS()
            self.model.write("MPmodel.ilp")
            
            
        self.obj=self.model.objval
        self.bound=self.model.objBound
        
        
        self.Opt_y={}
        self.Opt_beta1={}
        self.Opt_beta2={}
        
        for node in supply_nodes:
            self.Opt_y[node]={}
            for t in range(T):
                self.Opt_y[node][t]=node.y[t].X
                
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

        for node in supply_nodes :
            for t in range(T):
                self.model.addConstr((gu.quicksum(node.y[k] for k in range(t+1)) +node.In_cap >=node.MP_h[It][t]) , name="P_Cap"+str(It)+','+str(node.node_id)+','+str(t))
                self.model.addConstr((gu.quicksum(line.MP_z[It][t] for line in node.outgoing) <= node.MP_h[It][t]), name="Production"+str(It)+','+str(node.node_id)+','+str(t) )
        
        for node in port_nodes :
            for t in range(T):
                self.model.addConstr((gu.quicksum(line.MP_z[It][t] for line in node.outgoing) <= node.MP_v[It][t] ), name="Import"+str(It)+','+str(node.node_id)+','+str(t) )

        for node in demand_nodes :
            for t in range(T):
                self.model.addConstr((gu.quicksum(line.MP_z[It][t] for line in node.incoming) >= XI[node][t] ), name="Demad"+str(It)+','+str(node.node_id)+','+str(t) )

        self.model.addConstr(self.alpha+gu.quicksum((node.beta1[t]-node.beta2[t])*XI[node][t] for t in range(T) for node in demand_nodes) >= gu.quicksum(node.c_VP[t]*node.MP_h[It][t] for node in supply_nodes for t in range(T)) + gu.quicksum(node.c_I[t]*node.MP_v[It][t] for node in port_nodes for t in range(T) ) + gu.quicksum(-node.R[t]* XI[node][t] for node in demand_nodes for t in range(T)) + gu.quicksum(line.c_VT[t]*line.dist*line.MP_z[It][t] for line in Lines for t in range(T)), name='cut'+str(It) )
        
    def BendersCut(self,It,XI,SP_nu,SP_psi):
        self.model.addConstr(self.alpha+gu.quicksum((node.beta1[t]-node.beta2[t])*XI[node][t] for t in range(T) for node in demand_nodes) >= -gu.quicksum( gu.quicksum( node.y[k] for k in range(t+1))* SP_nu[node][t] for node in supply_nodes for t in range(T)) + gu.quicksum( XI[node][t]*(SP_psi[node][t]-node.R[t]) for node in demand_nodes for t in range(T)) ) 
    
    def UpdateObj(self):
        self.opening=gu.quicksum(node.c_FS[0]*node.x[0] for node in supply_nodes )+ gu.quicksum(node.c_FS[t]*(node.x[t]-node.x[t-1]) for node in supply_nodes for t in range(1,T))
        self.capacity=gu.quicksum(node.c_HS[t]*node.y[t] for node in supply_nodes for t in range(T))
        self.dual_beta1=gu.quicksum(node.mu[t]*node.beta1[t] +node.mu[t]*gu.quicksum(node.Lambda[supply]*node.phi1_variables[supply][t] for supply in supply_nodes) +node.epsilon[t]*node.beta1[t] for node in demand_nodes for t in range(T))
        self.dual_beta2=gu.quicksum(-node.mu[t]*node.beta2[t] -node.mu[t]*gu.quicksum(node.Lambda[supply]*node.phi2_variables[supply][t] for supply in supply_nodes) +node.epsilon[t]*node.beta2[t] for node in demand_nodes for t in range(T))
        self.total= self.opening+self.capacity + self.alpha +self.dual_beta1+ self.dual_beta2    
        self.model.setObjective( self.total ,gu.GRB.MINIMIZE)

class BudgetedSubProblem:

    def __init__(self, budget):
        self.budget = budget * sum((node.Support[t][-1] - node.Support[t][0])/2 for node in demand_nodes for t in range(T))
        print(self.budget)
        self.model= gu.Model()  
        self.model.setParam('OutputFlag', 0)
        self.model.Params.MIPGap = 0.0001 
        self.model.params.NonConvex = 2

        self.model.setParam('MIPFocus', 2)

        for t in range(T):
            for node in supply_nodes :
                node.AddSupplyVariables_SP(self.model,t)
                            
            for node in port_nodes :
                node.AddPortVariables_SP(self.model,t)
            
            for node in demand_nodes :
                node.AddDemandVariables_SP(self.model,t)
                node.AddDemandVariables_Budget(self.model,t)
                
            for line in Lines :
                line.AddVariables_SP(self.model,t)
            
            for node in supply_nodes :
                node.AddSupplyConstraints_SP(self.model,t)
            
            for node in port_nodes:
                node.AddPortConstraints_SP(self.model,t)
            
            for node in demand_nodes:
                node.AddDemandConstraints_SP(self.model,t)
                node.AddDemandConstraints_Budget(self.model,t)
        
        self.model.addConstr(gu.quicksum(node.w[t] for t in range(T) for node in demand_nodes) <= self.budget)
        
        self.productioncost= gu.quicksum(node.c_VP[t]*node.h[t] for node in supply_nodes for t in range(T))
        self.importcost= gu.quicksum(node.c_I[t]*node.v[t] for node in port_nodes  for t in range(T))
        self.transportationcost=gu.quicksum(line.c_VT[t]*line.dist*line.z[t] for line in Lines  for t in range(T))
        self.betasum=gu.quicksum(-(Opt_beta1[node][t]-Opt_beta2[node][t])*node.xi[t] for node in demand_nodes  for t in range(T))
        self.totalreward= gu.quicksum(-node.R[t]* node.xi[t] for node in demand_nodes  for t in range(T))
        
        self.total=self.productioncost+self.importcost+self.transportationcost+self.betasum+self.totalreward
        self.model.setObjective(self.total,gu.GRB.MAXIMIZE)
                                         
        self.model.update()

    def optimize(self,time_limit):
        self.model.setParam('TimeLimit', time_limit)
        self.start = timeit.default_timer()
        self.model.optimize()
        self.stop = timeit.default_timer()
        self.GurobiTime=self.stop-self.start
        self.model.write("BudgetSPmodel.lp")
        
        if self.model.status == 3 or self.model.status == 4 :
            self.model.computeIIS()
            self.model.write("BudgetSPmodel.ilp")
            
        
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
            self.opt_psi[k]={}
            self.opt_nu[k]={}
            for t in range(T):
                self.XI[k][t] = {}
                self.opt_psi[k][t]={}
                self.opt_nu[k][t]={}  
                for node in demand_nodes:
                    self.XI[k][t][node]=node.xi[t].xn
                    self.opt_psi[k][t][node]=node.psi[t].xn
                for node in supply_nodes:
                    self.opt_nu[k][t][node]=node.nu[t].xn

  
    def ChangeRHS(self,Opt_y):
        for t in range(T):
            for node in supply_nodes:
                self.model.getConstrByName("P_Cap"+str(node.node_id)+','+str(t)).rhs= -(sum(Opt_y[node][k] for k in range(t+1)) +node.In_cap)
                self.model.getConstrByName('CS_S1C'+str(node.node_id)+','+str(t)).rhs= -(sum(Opt_y[node][k] for k in range(t+1)) +node.In_cap)
    
    def ChangeObjective(self,Opt_beta1, Opt_beta2):
        for t in range(T):
            for node in demand_nodes:
                node.xi[t].Obj=-(Opt_beta1[node][t]-Opt_beta2[node][t])-node.R[t]

Mstart = timeit.default_timer()
MP=MasterProblem()

  
start = timeit.default_timer()

Opt_y, Opt_beta1, Opt_beta2 = MP.optimize(TIME_LIMIT)
Mstop = timeit.default_timer()

SP = BudgetedSubProblem(budget)

UB = float("inf")
LB = -float("inf")
AIt=0
It=0
perc5 = float("inf") ; perc2= float("inf") ; perc1=float("inf")
time5=float("inf") ;time2=float("inf") ; time1=float("inf") 
In_tol=0.50
End_tol=0.001

while (UB-LB)>=abs(0.001*UB) and timeit.default_timer()-start<=TIME_LIMIT :
    if It!=0:
        Opt_y, Opt_beta1, Opt_beta2=MP.optimize(TIME_LIMIT - (timeit.default_timer()-start))
        MP.Start()

    LB = MP.obj
    SPstart = timeit.default_timer()
    
    Opt_XI={}
    Opt_nu={}
    Opt_psi={}
    time_limit = TIME_LIMIT - (timeit.default_timer()-start)
    if time_limit < 0:
        break
    SP.ChangeObjective(Opt_beta1, Opt_beta2)
    SP.ChangeRHS(Opt_y)
    SP.optimize(time_limit) 
    SPstop = timeit.default_timer()
    SP_sum=SP.obj[0]
    if SP_sum+(MP.obj-MP.opt_alpha)<=UB:
            UB=SP_sum+(MP.obj-MP.opt_alpha)
            optimal=Opt_y
    xi_val = 0
    for node in demand_nodes:
        Opt_XI[node]={}
        Opt_psi[node]={}
        for t in range(T):
            Opt_XI[node][t]=SP.XI[0][t][node]
            Opt_psi[node][t]=SP.opt_psi[0][t][node]
            xi_val += Opt_XI[node][t]
    for node in supply_nodes:
        Opt_nu[node]={}
        for t in range(T):
            Opt_nu[node][t]=SP.opt_nu[0][t][node]
            
    if mode== Mode.BendersDecomposition:
        MP.BendersCut(It,Opt_XI,Opt_nu,Opt_psi)
    else:
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


List=[rand_seed,int(budget*100),supply_count, demand_count, port_count,T,LB,UB,round(100*(UB-LB)/abs(UB),4),It,stop-start,perc5,time5,perc2,time2,perc1,time1,optimal]
with open('computational_budget.csv', 'a') as f_object: 
    writer_object = writer(f_object)
    writer_object.writerow(List)
    f_object.close()
