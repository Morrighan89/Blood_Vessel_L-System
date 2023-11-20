
from pyvista.utilities.helpers import vtk_points
import vtk
import pyvista as pv
import numpy as np
import math
import copy
import os
from tqdm import tqdm
#from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
import vtkmodules.numpy_interface as np2vtk
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkPolyData,
    vtkPolyLine,
    vtkGenericCell
)
import gmsh
from VesselInterpreter import (
    ReadPolyData,
    ExtractLine,
    radiusArrayName,
    parallelTransportNormalsArrayName
)
from itertools import chain

def create_gmsh_mesh_carving(clFileName,ofile="gmsh.msh"):
    baseCl=ReadPolyData(clFileName)
    temp = os.path.splitext(ofile)
    var = (os.path.basename(temp[0]), temp[1])
    name_ofile=var[0]

    numberOfLines=baseCl.GetNumberOfCells()
    numberOfPoints=baseCl.GetNumberOfPoints()
    print(f" Vessel segments: {numberOfLines}, Nodes: {numberOfPoints}")
    options=['t1.geo','-tol', '1.e-14','-setnumber', 'Geometry.OCCFixDegenerated', '1','-setnumber', 'Geometry.OCCFixSmallEdges', '1','-setnumber', 'Geometry.OCCFixSmallFaces', '1']
    gmsh.initialize(argv=options)
    gmsh.model.add("DFG 3D")
   # print(gmsh.model.get  r('Geometry.Tolerance'))
    newdims=[]
    
    #lines=[692,693,694,695,696,697,698,935,936,937,938,939,940,941,942,943,944,945,946,947,948,949,950,951,952,953,954,955,956,957,958,959,960,961,962,963,964,965,966,967,968,969,970,971,972,973,974,975,976,977,978,979,980,981,982,983,984,985,986,987,988,989,990,991,992,993,994,995,996,997,998,999,1000,1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1034,1035,1126,1127,1128,1129,1130,1131,1132,1133,1134,1135,1136,1137,1138,1139,1140,1141,1142,1143,1144,1145,1146,1147,1148,1149,1150,1151,1152,1153,1154,1155,1156,1157,1158,1159,1160,1161,1162,1163,1221,1542,1543,1544,1545,1546,1547,1548,1549,1550,1551,1552,1556,1557,1558,1559,1660,1661,1662,1663,1664,1665,1755,2407,2408,2409,2410,2411,2412,2413,2414,2415,2416,2417,2418,2419,2420,2421,2422,2423,2424,2425,2457,2458,2459,2460,2461,2462,2463,2464,2465,2466,2467,2468,2469,2470,2471,2472,2473,2474,2475,2476,2477,2478,2479,2480,2481,2482,2483,2484,2604,2605,2606,2615,2647,2648,2649,2650,2651,2652,2653,2654,2655,3166,3167,3168,3169,3170,3171,3172,3173,3174,3175,3176,3177,3178,3179,3180,3181,3182,3183,3184,3185,3186,3187,3188,3189,3190,3191,3192,3193,3202,3223,3224,3225,3226,3227,3228,3229,3230,3231,3232,3233,3234,3235,3236,3373,3374,3379,3382,3383,3384,3385,3386,3397,3399,3400,3401,3402,3408,3435,3436,3437,3440,3441,4040,4041,4042,4043,4044,4045,4046,4047,4048,4049,4050,4051,4052,4053,4054,4055,4056,4057,4058,4059,4060,4061,4066,4087,4088,4089,4090,4091,4092,4093,4094,4095,4096,4097,4098,4099,4100,4101,4250,4251,4252,4253,4255,4259,4260,4270,4271,4292,4293,4294,4295,4297,4476,4477,4907,4908,4909,4910,4911,4912,4913,4914,4915,4916,4917,4918,4919,4920,4921,4922,4923,4924,4925,4926,4927,4928,4929,4930,4931,4932,4933,4934,4935,4936,4937,4938,4939,4940,4948,4969,4970,4971,4972,4973,4974,4975,4976,4977,4978,4979,4980,5143,5144,5151,5152,5168,5189,5190,5204,5205,5206,5207,5210,5405,5751,5752,5756,5869,5870,5871,5872,5873,5874,5875,5876,5877,5878,5879,5880,5881,5882,5883,5884,5885,5886,5887,5888,5889,5890,5904,5941,5942,5943,5944,5945,5946,5947,5948,5949,5950,5951,5952,5953,5954,5955,5956,5957,5958,5959,6135,6144,6145,6146,6147,6159,6160,6161,6162,6165,6166,6167,6173,6203,6204,6205,6206,6207,6208,6398,6862,6863,6864,6865,6866,6867,6868,6869,6870,6871,6872,6873,6874,6875,6876,6877,6878,6879,6880,6881,6882,6883,6884,6885,6886,6896,6924,6925,6926,6927,6928,6929,6930,6931,6932,6933,6934,6935,6936,6937,6938,6939,6940,6941,6942,6943,7155,7164,7165,7166,7167,7175,7176,7177,7178,7179,7180,7181,7182,7237,7238,7239,7240,7241,7242,7243,7244,7245,7246,7247,7250,7251,7252,7463,7464,7465,7965,7966,7967,7968,7969,7970,7971,7972,7973,7974,7975,7976,7977,7978,7979,7980,7981,7982,7983,7984,7985,7986,7987,7988,7989,7990,7999,8016,8017,8018,8019,8020,8021,8022,8023,8024,8025,8026,8027,8028,8029,8030,8031,8032,8033,8034,8035,8036,8037,8038,8039,8040,8300,8301,8302,8303,8304,8316,8317,8318,8319,8320,8332,8333,8334,8335,8340,8341,8342,8343,8396,8397,8398,8399,8400,8401,8610,8616,8617,9038,9197,9198,9199,9200,9201,9202,9203,9204,9205,9206,9207,9208,9209,9210,9211,9212,9213,9244,9245,9246,9247,9248,9249,9250,9251,9252,9253,9254,9255,9256,9257,9258,9259,9260,9261,9262,9263,9264,9265,9266,9555,9556,9557,9558,9559,9563,9564,9565,9566,9584,9585,9586,9587,9588,9589,9590,9591,9592,9613,9614,9615,9616,9617,9618,9619,9620,9621,9622,9626,9627,9681,9682,9683,9684,9685,9686,9687,9688,9689,9690,9691,9692,9693,9694,9695,9696,9697,9917,9918,9921,9922,9923,10579,10580,10581,10582,10583,10584,10585,10586,10587,10588,10589,10590,10591,10592,10593,10594,10595,10620,10621,10622,10623,10624,10625,10626,10627,10628,10629,10630,10631,10632,10633,10634,10635,10636,10637,10638,10639,10640,10641,10642,10941,10955,10956,10957,10965,10966,10967,10983,10984,10985,10986,10987,10988,10989,11000,11001,11063,11064,11065,11066,11067,11068,11069,11070,11071,11072,11073,11350,11356,11357,11871,12042,12043,12044,12045,12046,12047,12048,12049,12050,12051,12052,12053,12054,12055,12085,12086,12087,12088,12089,12090,12091,12092,12093,12094,12095,12096,12097,12098,12099,12446,12447,12468,12473,12474,12475,12493,12494,12495,12496,12497,12509,12560,12561,12562,12563,12564,12565,12566,12567,12568,12569,12570,12571,12572,12573,12777,12778,12779,12780,12781,12784,13423,13424,13425,13426,13427,13428,13429,13430,13431,13432,13433,13434,13463,13464,13465,13466,13467,13468,13469,13470,13471,13472,13473,13474,13475,13476,13477,13478,13479,13480,13481,13482,13483,13830,13850,13853,13854,13872,13873,13874,13875,13876,13877,13878,13879,13880,13881,13882,13892,13951,13952,13953,13954,13955,13956,13957,13958,13959,13960,13972,13973,13974,13975,14195,14196,14766,14767,14768,14769,14770,14771,14772,14801,14802,14803,14804,14805,14806,14807,14808,14809,14810,14811,14812,14813,14814,15164,15165,15166,15167,15171,15172,15201,15202,15203,15230,15231,15232,15233,15234,15235,15236,15237,15254,15255,15256,15257,15318,15319,15320,15321,15548,15549,15550,15552,16161,16162,16163,16164,16192,16193,16194,16195,16196,16197,16198,16199,16200,16201,16202,16203,16204,16205,16206,16207,16208,16209,16210,16211,16212,16213,16214,16552,16553,16554,16556,16557,16558,16559,16576,16587,16609,16610,16611,16612,16613,16614,16615,16616,16628,16629,16630,16631,16632,16633,16634,16635,16636,16637,16716,16717,16718,16719,16720,16721,16726,16912,17537,17538,17560,17561,17562,17563,17564,17565,17566,17567,17568,17569,17570,17571,17572,17573,17574,17908,17909,17910,17911,17912,17913,17938,17951,17952,17984,17985,17986,17987,17988,17989,17990,17991,17997,17998,17999,18000,18001,18002,18003,18004,18005,18006,18007,18008,18009,18010,18011,18059,18060,18061,18786,18787,18808,18809,18810,18811,18812,18813,18814,18815,18816,18817,18818,18819,18820,18821,19150,19151,19168,19169,19175,19206,19214,19215,19216,19217,19218,19219,19287,19288,19289,19290,19291,19292,19462,19463,19967,19968,19969,19970,19987,19988,19989,19990,19991,19992,19993,19994,19995,19996,19997,19998,19999,20000,20001,20002,20003,20004,20005,20373,20374,20375,20426,20427,20428,20429,20430,20431,20432,20433,20434,20435,21024,21025,21052,21053,21054,21055,21056,21057,21058,21059,21060,21061,21062,21063,21064,21065,21066,21067,21068,21069,21407,21408,21409,21410,21411,21474,21475,21476,21477,21478,21479,21480,21481,22055,22056,22057,22058,22059,22060,22061,22062,22063,22383,22384,22385,22386,22387,22388,22389,22453,22454,22455,22456,22457,22936,22954,22955,22956,23249,23250,23251,23252,23253,23254,23255,23309,23310,23311,23312,23313,23314,23315,23316,23730,23731,23732,23733,23734,23735,23736,23737,23738,24017,24018,24019,24020,24072,24073,24074,24075,24076,24448,24449,24450,24730,25028,25029,25030,25219,25220,25257,25258,25259,25260,25261,25262,25531,25532,25533,25691,25692,25693,25694,25695,25696,25697,25868,25869,25989,26024,26025,26026,26027,26291,26292,26293,26294,26295,26296,26297,26298,26299,26300,26301,26470,26471,26472,26473,26593,26594,26595]
    #cx,cy,cz=0.8734448800133734, -2.081709684995058, 1.6177133644915518
    #dx,dy,dz=1, 1, 1
    #cx,cy,cz=0.2816580241399551, -1.9043948378959654, -1.3874804023695808
    #dx,dy,dz=3.195877166439165,2.9768761114761113,2.5549797529430096
    #cx,cy,cz=0.22723123693223046,-1.117517697735693,-1.5175201865734897  # villi14
    #dx,dy,dz=2.0870246230892904,2.8587875728566874,3.0326460420042394    # villi14
    cx,cy,cz=0.4897231334111418,-0.4990486324439569,-1.4965857674994025# villi15
    dx,dy,dz=1.9164413473292243,2.5604400762919113,2.7007688261857514# villi15
    
    box=gmsh.model.occ.addBox(cx,cy,cz, dx,dy,dz)
    volumes=gmsh.model.occ.getEntities(dim=3)
    count = 0
    for i in tqdm(range(0,numberOfLines)):#numberOfLineschain(range(0,136),range(138,numberOfLines))
        line=ExtractLine(i,baseCl)
        linePtsNumber=line.GetNumberOfPoints()
        newtrunks=[]
        volumes=gmsh.model.occ.getEntities(dim=3)
        try:
        #volumes=gmsh.model.occ.getEntities(dim=3) 
            for j in range(0,linePtsNumber-1): #linePtsNumber-1
                startPoint=line.GetPoint(j)
                startingRadius=line.GetPointData().GetArray(radiusArrayName).GetTuple1(j)
                endingRadius=line.GetPointData().GetArray(radiusArrayName).GetTuple1(j+1)
                direction=line.GetPointData().GetArray(parallelTransportNormalsArrayName).GetTuple3(j)
                if startingRadius==endingRadius:
                    cone = gmsh.model.occ.addCylinder(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius)
                else:
                    cone = gmsh.model.occ.addCone(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius,endingRadius)
                    newtrunks.append((3,cone))
                if j<linePtsNumber-2:
                    ball = gmsh.model.occ.addSphere(startPoint[0]+direction[0],startPoint[1]+direction[1],startPoint[2]+direction[2],endingRadius)
                    newtrunks.append((3,ball))
                #if j==0:
                #    ball = gmsh.model.occ.addSphere(startPoint[0],startPoint[1],startPoint[2],startingRadius*1.01)
                #    newtrunks.append((3,ball))
            #volumes=gmsh.model.occ.getEntities(dim=3)
        
            #pippo=gmsh.model.occ.fuse([newtrunks[0]], newtrunks)
            gmsh.model.occ.cut(volumes,newtrunks,removeTool= True)
            gmsh.model.occ.synchronize()
            #gmsh.model.occ.cut(volumes,pippo[0],removeObject= True, removeTool= True) 
        except:
            print(f'Problem with polyline {i}')
        volumes=gmsh.model.occ.getEntities(dim=3)
        gmsh.model.occ.synchronize()
        count=count+1
        if ((count % 1) == 0 or i==1):
            gmsh.write(f'temp\model_{name_ofile}_{i}.brep')
    #cx,cy,cz=0.8734448800133734, -2.081709684995058, 1.6177133644915518
    #cx,cy,cz=0.12864125650976543,-2.0758902928179417,-1.8674509158982868
    #dx,dy,dz=2.580985600142573,3.589531105579775,3.400352085492524
    #dx,dy,dz=1, 1, 1
    #box=gmsh.model.occ.addBox(cx,cy,cz, dx,dy,dz)
    #gmsh.model.occ.cut([(3,box)],volumes, removeTool= True) 
    gmsh.model.occ.synchronize()
    gmsh.write(f'temp\model_{name_ofile}_fine.brep')
    print(f'Begin labeling')
    surfaces = gmsh.model.occ.getEntities(dim=2)
    root_marker, leaves_marker, wall_marker, bottom_side_marker, top_side_marker, lateral_sides_marker = 203, 205, 207, 209, 211, 213
    walls = []
    leaves = []
    laterals=[]
    comtop =gmsh.model.occ.getCenterOfMass(2, 1)
    for surface in tqdm(surfaces):
        com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
        surface_type= gmsh.model.getType(surface[0], surface[1])
        if surface_type == 'Plane':
            if np.allclose(com, comtop):
                gmsh.model.addPhysicalGroup(surface[0], [surface[1]], top_side_marker)
                gmsh.model.setPhysicalName(surface[0], top_side_marker, "topSide")
            elif np.allclose(com, [cx+dx, cy+dy/2, cz+dz/2]):
                gmsh.model.addPhysicalGroup(surface[0], [surface[1]], bottom_side_marker)
                gmsh.model.setPhysicalName(surface[0], bottom_side_marker, "bottomSide")
            elif np.allclose(com, [0,0,0]):
                gmsh.model.addPhysicalGroup(surface[0], [surface[1]], root_marker)
                gmsh.model.setPhysicalName(surface[0], root_marker, "root")
            elif np.isclose(com[2], cz) or np.isclose(com[1], cy+dy) or np.isclose(com[2], cz+dz) or np.isclose(com[1],cy):
                laterals.append(surface[1])
            else:
                leaves.append(surface[1])
        else:
            walls.append(surface[1])
    gmsh.model.addPhysicalGroup(2, laterals, lateral_sides_marker)
    gmsh.model.setPhysicalName(2, lateral_sides_marker, "lateral")
    gmsh.model.addPhysicalGroup(2, leaves, leaves_marker)
    gmsh.model.setPhysicalName(2, leaves_marker, "leaves")
    gmsh.model.addPhysicalGroup(2, walls, wall_marker)
    gmsh.model.setPhysicalName(2, wall_marker, "Walls")
    gmsh.model.occ.synchronize()
    gmsh.write(f'model_{name_ofile}.brep')
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0005)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.18)
    gmsh.model.mesh.generate(3)
    gmsh.write(f'model_{name_ofile}.msh')

def create_gmsh_mesh_adding(clFileName,folder="temp",ofile="gmsh.msh",scalingfactor=1):
    baseCl=ReadPolyData(clFileName)
    temp = os.path.splitext(ofile)
    var = (os.path.basename(temp[0]), temp[1])
    name_ofile=var[0]
    numberOfLines=baseCl.GetNumberOfCells()
    numberOfPoints=baseCl.GetNumberOfPoints()
    print(f" Vessel segments: {numberOfLines}, Nodes: {numberOfPoints}")
    options=['t1.geo','-tol', '1.e-9','-setnumber', 'Geometry.OCCFixDegenerated', '1','-setnumber', 'Geometry.OCCFixSmallEdges', '1','-setnumber', 'Geometry.OCCFixSmallFaces', '1']
    gmsh.initialize(argv=options)
    gmsh.initialize()
    gmsh.model.add("DFG 3D")
    
    newdims=[]
    line=ExtractLine(0,baseCl)
    startPoint=line.GetPoint(0)
    startingRadius=scalingfactor*line.GetPointData().GetArray(radiusArrayName).GetTuple1(0)
    endingRadius=scalingfactor*line.GetPointData().GetArray(radiusArrayName).GetTuple1(1)    
    direction=line.GetPointData().GetArray(parallelTransportNormalsArrayName).GetTuple3(0)
    if startingRadius==endingRadius:
        trunk = gmsh.model.occ.addCylinder(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius)
    else:
        trunk = gmsh.model.occ.addCone(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius, endingRadius)
    ball = gmsh.model.occ.addSphere(startPoint[0]+direction[0],startPoint[1]+direction[1],startPoint[2]+direction[2],endingRadius)
    for i in range(14, numberOfLines+1):#chain(range(381,401),range(401,numberOfLines)):#range(0, numberOfLines):#chain(range(0,136),range(138,161),range(181,numberOfLines)):#chain(range(0,136),range(138,numberOfLines))
        line=ExtractLine(i,baseCl)
        linePtsNumber=line.GetNumberOfPoints()
        newtrunks=[]

        volumes=gmsh.model.occ.getEntities(dim=3) 
        for j in range(1,linePtsNumber-1): #linePtsNumber-1
            startPoint=line.GetPoint(j)
            startingRadius=scalingfactor*line.GetPointData().GetArray(radiusArrayName).GetTuple1(j)
            endingRadius=scalingfactor*line.GetPointData().GetArray(radiusArrayName).GetTuple1(j+1)
            direction=line.GetPointData().GetArray(parallelTransportNormalsArrayName).GetTuple3(j)
            if startingRadius==endingRadius:
                cone = gmsh.model.occ.addCylinder(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius)
            else:
                cone = gmsh.model.occ.addCone(startPoint[0],startPoint[1],startPoint[2],direction[0], direction[1],direction[2], startingRadius,endingRadius)
                newtrunks.append((3,cone))
            if j<linePtsNumber-2:
                ball = gmsh.model.occ.addSphere(startPoint[0]+direction[0],startPoint[1]+direction[1],startPoint[2]+direction[2],endingRadius*1.01)
                newtrunks.append((3,ball))
            #if j==0:
            #    ball = gmsh.model.occ.addSphere(startPoint[0],startPoint[1],startPoint[2],startingRadius*1.01)
            #    newtrunks.append((3,ball))
        volumes=gmsh.model.occ.getEntities(dim=3)
        gmsh.model.occ.fuse(newtrunks,volumes,removeObject= True, removeTool= True)
        #gmsh.model.occ.fuse([(3,trunk)], newtrunks,removeObject= False, removeTool= True)
        volumes=gmsh.model.occ.getEntities(dim=3)
        print(volumes,i)
        gmsh.model.occ.synchronize()

        if ((i % 1) == 0  or i==1):
            gmsh.write(f'{folder}\{name_ofile}_{i}.step')
    #cx,cy,cz=0.2816580241399551, -1.9043948378959654, -1.3874804023695808
    #dx,dy,dz=3.195877166439165,2.9768761114761113,2.5549797529430096
    #cx,cy,cz=0.22723123693223046,-1.117517697735693,-1.5175201865734897  # villi14
    #dx,dy,dz=2.0870246230892904,2.8587875728566874,3.0326460420042394    # villi14
    #cx,cy,cz=0.11962770121213051,-0.4990486324439569,-1.4965857674994025# villi15
    #dx,dy,dz=2.2865367795282356,2.5604400762919113,2.7007688261857514# villi15
    #cx,cy,cz=0,-1.2402858552276663,-1.3034999407295926# villi16
    #dx,dy,dz=2.2,2.4,3# villi16
    cx,cy,cz=0.2, -1.75, -1.75 # villi19
    dx,dy,dz=2.2, 3.5, 3.5 # villi19
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0005)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.18)
    gmsh.model.mesh.generate(2)
    gmsh.write(f'{folder}\model_{name_ofile}_2.vtk')
    box=gmsh.model.occ.addBox(cx,cy,cz, dx,dy,dz)
    
    gmsh.model.occ.cut([(3,box)],volumes, removeTool= True) 
    gmsh.model.occ.synchronize()
    gmsh.write(f'{folder}\model_{name_ofile}.step')
    print(f'Begin labeling')
    surfaces = gmsh.model.occ.getEntities(dim=2)
    root_marker, leaves_marker, wall_marker, bottom_side_marker, top_side_marker, lateral_sides_marker = 203, 205, 207, 209, 211, 213
    walls = []
    leaves = []
    laterals=[]
    volumes = gmsh.model.occ.getEntities(dim=3)
    gmsh.model.addPhysicalGroup(3, [2], 300)
    gmsh.model.setPhysicalName(3, 300, "volume")
    print(volumes)
    comtop =gmsh.model.occ.getCenterOfMass(2, 1)
    for surface in tqdm(surfaces):
        com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
        surface_type= gmsh.model.getType(surface[0], surface[1])
        if surface_type == 'Plane':
            if np.allclose(com, comtop):
                gmsh.model.addPhysicalGroup(surface[0], [surface[1]], top_side_marker)
                gmsh.model.setPhysicalName(surface[0], top_side_marker, "topSide")
            elif np.allclose(com, [cx+dx, cy+dy/2, cz+dz/2]):
                gmsh.model.addPhysicalGroup(surface[0], [surface[1]], bottom_side_marker)
                gmsh.model.setPhysicalName(surface[0], bottom_side_marker, "bottomSide")
            elif np.allclose(com, [0,0,0]):
                gmsh.model.addPhysicalGroup(surface[0], [surface[1]], root_marker)
                gmsh.model.setPhysicalName(surface[0], root_marker, "root")
            elif np.isclose(com[2], cz) or np.isclose(com[1], cy+dy) or np.isclose(com[2], cz+dz) or np.isclose(com[1],cy):
                laterals.append(surface[1])
            else:
                leaves.append(surface[1])
        else:
            walls.append(surface[1])
    gmsh.model.addPhysicalGroup(2, laterals, lateral_sides_marker)
    gmsh.model.setPhysicalName(2, lateral_sides_marker, "lateral")
    gmsh.model.addPhysicalGroup(2, leaves, leaves_marker)
    gmsh.model.setPhysicalName(2, leaves_marker, "leaves")
    gmsh.model.addPhysicalGroup(2, walls, wall_marker)
    gmsh.model.setPhysicalName(2, wall_marker, "Walls")
    gmsh.model.occ.synchronize()
    gmsh.write(f'{folder}\model_{name_ofile}_2.step')
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0005)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.18)
    gmsh.model.mesh.generate(3)
    gmsh.write(f'{folder}\model_{name_ofile}_3.msh')
    gmsh.write(f'{folder}\model_{name_ofile}_3.vtk')
    gmsh.finalize()
    
def boolean_example():
    gmsh.initialize()

    gmsh.model.add("boolean")

    # from http://en.wikipedia.org/wiki/Constructive_solid_geometry

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.4)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.4)

    R = 1.4
    Rs = R * .7
    Rt = R * 1.25

    gmsh.model.occ.addBox(-R, -R, -R, 2 * R, 2 * R, 2 * R, 1)
    gmsh.model.occ.addSphere(0, 0, 0, Rt, 2)
    gmsh.model.occ.intersect([(3, 1)], [(3, 2)], 3)
    gmsh.model.occ.addCylinder(-2 * R, 0, 0, 4 * R, 0, 0, Rs, 4)
    gmsh.model.occ.addCylinder(0, -2 * R, 0, 0, 4 * R, 0, Rs, 5)
    gmsh.model.occ.addCylinder(0, 0, -2 * R, 0, 0, 4 * R, Rs, 6)
    gmsh.model.occ.fuse([(3, 4), (3, 5)], [(3, 6)], 7)
    gmsh.model.occ.cut([(3, 3)], [(3, 7)], 8)

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.generate(3)
    #gmsh.model.mesh.refine()
    #gmsh.model.mesh.setOrder(2)
    #gmsh.model.mesh.partition(4)

    gmsh.write("boolean.msh")
    gmsh.write("boolean.brep")
    gmsh.finalize()

def example_spline_extrude():
    gmsh.initialize()

    gmsh.model.add("extrude spline")
    nturns = 2 # tested ok up to 100

    npts = 100 * nturns
    r = 1.
    rd = 0.1
    h = 1. * nturns

    for i in range(npts):
      theta = i * 2. * math.pi * nturns / npts
      gmsh.model.occ.addPoint(r * math.cos(theta), r * math.sin(theta),
                              i * h / npts, i+1)

    gmsh.model.occ.addSpline(range(1, npts), 1)
    gmsh.model.occ.addWire([1], 1)

    gmsh.model.occ.addDisk(1,0,0, rd, rd, 1)

    gmsh.model.occ.addRectangle(1+2*rd,-rd,0, 2*rd,2*rd, 2, rd/5)
    gmsh.model.occ.rotate([(2, 1), (2, 2)], 0, 0, 0, 1, 0, 0, math.pi/2)

    #gmsh.model.occ.addPipe([(2, 1), (2, 2)], 1, 'DiscreteTrihedron')
    gmsh.model.occ.addPipe([(2, 1), (2, 2)], 1, 'Frenet')

    gmsh.model.occ.remove([(2, 1), (2, 2), (1, 1)])

    gmsh.model.occ.synchronize()

    gmsh.option.setNumber('Mesh.MeshSizeMin', 0.1)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 0.1)
    gmsh.option.setNumber('Geometry.NumSubEdges', npts) # nicer display of curves
    gmsh.model.mesh.generate(3)
    gmsh.write("spline.msh")
    gmsh.write("spline.brep")
    gmsh.finalize()

def main():# used for debug purpose of function in this file
    create_gmsh_mesh_adding("vtkVilli16.vtp",folder="step",ofile='villi16.msh',scalingfactor=1)
    #create_gmsh_mesh_carving("vtkVilli10Trunc.vtp")
    #create_gmsh_mesh_carving("vtkVilli15.vtp",'villi15.msh')
    #boolean_example()
    #example_spline_extrude()
if __name__=='__main__':
    
    main()