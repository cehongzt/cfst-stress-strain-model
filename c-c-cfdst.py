#!/user/bin/python
#-*-coding:UTF-8-*-
from abaqus import *
from abaqusConstants import *
from caeModules import *
import math
logical=True
while logical:
	parameter=(('Model Name:','Model-Circular-CFDST-1'),('Diameter_OUT(mm):','300'),('Thickness_OUT(mm):','8'),('Diameter_IN(mm):','200'),('Thickness_IN(mm):','8'),
		('Length(mm):','800'),('fy_OUT(N/mm):','345'),('fy_IN(N/mm):','345'),('fcu(N/mm):','30'),('Displacement(mm):','30'),
		('Mesh number(mm):','8'),('Job Name:','Job-Circular-CFDST-01'))
	Modelname,diameter_out,thickness_out,diameter_in,thickness_in,length,fyo,fyi,fcu,displacement,meshsize,jobname=getInputs(fields=parameter,
		label='Please provide the following information:',dialogTitle='Model Parameter')
	listtry=[Modelname,diameter_out,thickness_out,diameter_in,thickness_in,length,fyo,fyi,fcu,displacement,meshsize,jobname]

	if listtry==[None,None,None,None,None,None,None,None,None,None,None,None]:
		logical=False
	else:
		logical=False
        viewportshow = session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=172.5, height=151.0)
        session.viewports['Viewport: 1'].makeCurrent()
        session.viewports['Viewport: 1'].maximize()
        diameter_out=float(diameter_out)
        thickness_out=float(thickness_out)
        diameter_in=float(diameter_in)
        thickness_in=float(thickness_in)
        length=float(length)
        fyo=float(fyo)
        fyi=float(fyi)
        fcu=float(fcu)

        displacement=float(displacement)
        meshsize=int(meshsize)
        sheetsize=diameter_out+50
##������part
        myModel=mdb.Model(name=Modelname)
        s= myModel.ConstrainedSketch(name='A',sheetSize=sheetsize)
        s.setPrimaryObject(option=STANDALONE)
        s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=((diameter_out/2-thickness_out), 0.0))
        s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=((diameter_in / 2 ), 0.0))
        p = myModel.Part(name='Part-concrete', dimensionality=THREE_D,type=DEFORMABLE_BODY)
        p = myModel.parts['Part-concrete']
        p.BaseSolidExtrude(sketch=s, depth=length)
##���part
        s = myModel.ConstrainedSketch(name='A',sheetSize=sheetsize)
        s.setPrimaryObject(option=STANDALONE)
        s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(diameter_out/2, 0.0))
        p = myModel.Part(name='Part-rigid', dimensionality=THREE_D,type=DEFORMABLE_BODY)
        p = myModel.parts['Part-rigid']
        p.BaseSolidExtrude(sketch=s, depth=20.0)

        s = myModel.ConstrainedSketch(name='A',sheetSize=sheetsize)
        s.setPrimaryObject(option=STANDALONE)
        s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=((diameter_out/2-thickness_out), 0.0))
        p = myModel.Part(name='Part-tube-out', dimensionality=THREE_D,type=DEFORMABLE_BODY)
        p = myModel.parts['Part-tube-out']
        p.BaseShellExtrude(sketch=s, depth=length)

        s = myModel.ConstrainedSketch(name='A', sheetSize=sheetsize)
        s.setPrimaryObject(option=STANDALONE)
        s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=((diameter_in / 2 - thickness_in), 0.0))
        p = myModel.Part(name='Part-tube-in', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p = myModel.parts['Part-tube-in']
        p.BaseShellExtrude(sketch=s, depth=length)
        s.unsetPrimaryObject()
        del myModel.sketches['A']
        myModel.Material(name='Material-rigid')
        myModel.materials['Material-rigid'].Density(table=((7.8e-09,),))
        myModel.materials['Material-rigid'].Elastic(table=((1000000000000.0, 1e-06),))

        Es=206000

        epsilone = 0.8 * fyo / Es
        epsilone1 = 1.5 * epsilone
        epsilone2 = 10 * epsilone1
        epsilone3 = 100 * epsilone1
        A = 0.2 * fyo / (epsilone1 - epsilone) / (epsilone1 - epsilone)
        B = 2 * A * epsilone1
        C = 0.8 * fyo + A * epsilone * epsilone - B * epsilone
        list1 = []
        tuple1 = ()



        def drange(start, stop, step):
            r = start
            while r < (stop - step):
                yield r
                r += step



        for epsilons in drange(epsilone * 100000, epsilone1 * 100000, 5):
            sigma = -A * epsilons * epsilons / 10000000000 + B * epsilons / 100000 + C
            epsilon = (epsilons / 100000) - epsilone
            sigma = float('%.2f' % sigma)
            epsilon = float('%.6f' % epsilon)
            list1.append((sigma, epsilon))

        for epsilons in drange(epsilone1 * 10000, epsilone2 * 10000, 5):
            sigma = fyo
            epsilon = (epsilons / 10000) - epsilone
            sigma = float('%.2f' % sigma)
            epsilon = float('%.6f' % epsilon)
            list1.append((sigma, epsilon))

        for epsilons in drange(epsilone2 * 1000, epsilone3 * 1000, 5):
            sigma = fyo * (1 + 0.6 * ((epsilons - epsilone2 * 1000) / (epsilone3 * 1000 - epsilone2 * 1000)))
            epsilon = (epsilons / 1000) - epsilone
            sigma = float('%.2f' % sigma)
            epsilon = float('%.6f' % epsilon)
            list1.append((sigma, epsilon))

        for epsilons in drange(epsilone3 * 1000, epsilone3 * 2000, 50):
            sigma = 1.6 * fyo
            epsilon = (epsilons / 1000) - epsilone
            sigma = float('%.2f' % sigma)
            epsilon = float('%.6f' % epsilon)
            list1.append((sigma, epsilon))
        tuple1=tuple(list1)

        steel_out = myModel.Material(name='Material-steel-out')
        steel_out.Density(table=((7.8e-9,),))
        steel_out.Elastic(table=((206000.0, 0.3),))
        steel_out.Plastic(table=tuple1)

        epsilone_in = 0.8 * fyi / Es
        epsilone_in1 = 1.5 * epsilone_in
        epsilone_in2 = 10 * epsilone_in1
        epsilone_in3 = 100 * epsilone_in1
        A1 = 0.2 * fyi / (epsilone_in1 - epsilone_in) / (epsilone_in1 - epsilone_in)
        B1 = 2 * A1 * epsilone_in1
        C1 = 0.8 * fyi + A1 * epsilone_in * epsilone_in - B1 * epsilone_in
        list3 = []
        tuple3= ()






        for epsilonsin in drange(epsilone_in * 100000, epsilone_in1 * 100000, 5):
            sigma_in = -A1 * epsilonsin * epsilonsin / 10000000000 + B1 * epsilonsin / 100000 + C1
            epsilon_in  = (epsilonsin / 100000) - epsilone_in
            sigma_in = float('%.2f' % sigma_in)
            epsilon_in  = float('%.6f' % epsilon_in )
            list3.append((sigma_in, epsilon_in ))

        for epsilonsin in drange(epsilone_in1 * 10000, epsilone_in2 * 10000, 5):
            sigma_in = fyi
            epsilon_in  = (epsilonsin / 10000) - epsilone_in
            sigma_in = float('%.2f' % sigma_in)
            epsilon_in  = float('%.6f' % epsilon_in )
            list3.append((sigma_in, epsilon_in ))

        for epsilonsin in drange(epsilone_in2 * 1000, epsilone_in3 * 1000, 5):
            sigma_in = fyi * (1 + 0.6 * ((epsilonsin - epsilone_in2 * 1000) / (epsilone_in3 * 1000 - epsilone_in2 * 1000)))
            epsilon_in  = (epsilonsin / 1000) - epsilone_in
            sigma_in = float('%.2f' % sigma_in)
            epsilon_in  = float('%.6f' % epsilon_in )
            list3.append((sigma_in, epsilon_in ))

        for epsilonsin in drange(epsilone_in3 * 1000, epsilone_in3 * 2000, 50):
            sigma_in = 1.6 * fyi
            epsilon_in  = (epsilonsin / 1000) - epsilone_in
            sigma_in = float('%.2f' % sigma_in)
            epsilon_in  = float('%.6f' % epsilon_in )
            list3.append((sigma_in, epsilon_in ))
        tuple3=tuple(list3)
        steel_in = myModel.Material(name='Material-steel-in')
        steel_in.Density(table=((7.8e-9,),))
        steel_in.Elastic(table=((206000.0, 0.3),))
        steel_in.Plastic(table=tuple3)

        if fcu <= 50:
            alpha1 = 0.76
        elif fcu >= 80:
            alpha1 = 0.82
        else:
            alpha1 = 0.76 + 0.002 * (fcu - 50)

        if fcu <= 40:
            alpha2 = 1
        elif fcu >= 80:
            alpha2 = 0.87
        else:
            alpha2 = 1 + 0.00325 * (40 - fcu)
        fc = 0.79 * fcu
        fck = 0.88 * alpha1 * alpha2 * fcu

        Ec = 4700 * sqrt(fc)
        miu = 0.3

        Aso = 0.25 * pi * ((diameter_out ** 2) - ((diameter_out - 2 * thickness_out) ** 2))
        Ace = 0.25 * pi * ((diameter_out - 2 * thickness_out) ** 2)
        ksi = Aso * fyo / Ace / fck
        #
        list2 = []
        tuple2 = ()
        eta = 2  #
        sigma0 = fc
        epsilonc = (1300 + 12.5 * fc) * 10 ** (-6)
        epsilon0 = epsilonc + 800 * (ksi ** 0.2) * (10 ** (-6))
        beta0 = 0.5 * (((2.36 * (10 ** (-5))) ** (0.25 + ((ksi - 0.5) ** 7))) * (fc ** 0.5))
        if beta0 < 0.12:
            beta0 = 0.12
        concretestart = 0.3 * fc / Ec
        for epsiloncc in drange(concretestart, concretestart + 25 * 0.0002, 0.0002):
            x = epsiloncc / epsilon0
            if x <= 1:
                y = 2 * x - x ** 2
            if x > 1:
                y = x / ((beta0 * ((x - 1) ** eta)) + x)
            sigmac = y * sigma0
            epsilonccc = epsiloncc - concretestart
            sigmac = float('%.2f' % sigmac)
            epsilonccc = float('%.6f' % epsilonccc)
            list2.append((sigmac, epsilonccc))
        for epsiloncc in drange(concretestart + 25 * 0.0002, 0.1, 0.0005):
            x = epsiloncc / epsilon0
            if x <= 1:
                y = 2 * x - x ** 2
            if x > 1:
                y = x / ((beta0 * ((x - 1) ** eta)) + x)
            sigmac = y * sigma0
            epsilonccc = epsiloncc - concretestart
            sigmac = float('%.2f' % sigmac)
            epsilonccc = float('%.6f' % epsilonccc)
            list2.append((sigmac, epsilonccc))
        tuple2 = tuple(list2)
        # ??????????????????
        sigmat0 = 0.26 * ((1.25 * fc) ** 0.666667)
        sigmat0 = float('%.3f' % sigmat0)
        if fcu == 20:
            gfi = 0.04
        elif fcu == 40:
            gfi = 0.12
        else:
            gfi = 0.004 * fcu - 0.04
        EC = float('%.2f' % Ec)

        concrete = myModel.Material(name='Material-concrete')
        concrete.Density(table=((2.4e-9,),))
        concrete.Elastic(table=((EC, 0.2),))
        concrete.ConcreteDamagedPlasticity(table=((30.0, 0.1, 1.16, 0.66667, 1e-05),))
        concrete.concreteDamagedPlasticity.ConcreteCompressionHardening(table=(tuple2))
        concrete.concreteDamagedPlasticity.ConcreteTensionStiffening(table=((sigmat0, gfi),), type=GFI)
        myModel.HomogeneousSolidSection(name='Section-concrete', material='Material-concrete', thickness=None)
        myModel.HomogeneousSolidSection(name='Section-rigid', material='Material-rigid', thickness=None)
        myModel.HomogeneousShellSection(name='Section-tube-in',
                                        preIntegrate=OFF, material='Material-steel-in', thicknessType=UNIFORM,
                                        thickness=thickness_in, thicknessField='', idealization=NO_IDEALIZATION,
                                        poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT,
                                        useDensity=OFF, integrationRule=SIMPSON, numIntPts=9)
        myModel.HomogeneousShellSection(name='Section-tube-out',
                                        preIntegrate=OFF, material='Material-steel-out', thicknessType=UNIFORM,
                                        thickness=thickness_out, thicknessField='', idealization=NO_IDEALIZATION,
                                        poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT,
                                        useDensity=OFF, integrationRule=SIMPSON, numIntPts=9)
        p = myModel.parts['Part-concrete']

        c = p.cells
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        region = regionToolset.Region(cells=cells)
        p.SectionAssignment(region=region, sectionName='Section-concrete', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)
        p = myModel.parts['Part-rigid']
        c = p.cells
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        region = regionToolset.Region(cells=cells)
        p.SectionAssignment(region=region, sectionName='Section-rigid', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)
        p = myModel.parts['Part-tube-out']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region = regionToolset.Region(faces=faces)
        p.SectionAssignment(region=region, sectionName='Section-tube-out', offset=0.0,
            offsetType=BOTTOM_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)
        p = myModel.parts['Part-tube-in']
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        region = regionToolset.Region(faces=faces)
        p.SectionAssignment(region=region, sectionName='Section-tube-in', offset=0.0,
            offsetType=BOTTOM_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)

        # Assembly???
        a = myModel.rootAssembly
        session.viewports['Viewport: 1'].setValues(displayedObject=a)
        a.DatumCsysByDefault(CARTESIAN)
        p = myModel.parts['Part-concrete']
        a.Instance(name='Part-concrete-1', part=p, dependent=OFF)
        p = myModel.parts['Part-rigid']
        a.Instance(name='Part-rigid-1', part=p, dependent=OFF)
        a.Instance(name='Part-rigid-2', part=p, dependent=OFF)
        p = myModel.parts['Part-tube-out']
        a.Instance(name='Part-tube-1-out', part=p, dependent=OFF)
        p = myModel.parts['Part-tube-in']
        a.Instance(name='Part-tube-2-in', part=p, dependent=OFF)


        ##????????????
        a.translate(instanceList=('Part-rigid-1',), vector=(0.0, 0.0, -20.0))  # 30???????
        a.translate(instanceList=('Part-rigid-2',), vector=(0.0, 0.0, length))  # 800?????????
        #step
        myModel.StaticStep(name='Step-1', previous='Initial',
                           maxNumInc=1000000, initialInc=0.0001, minInc=1e-06, maxInc=0.1, nlgeom=ON)
        ##????set
        f1 = a.instances['Part-rigid-1'].faces
        faces1 = f1.getSequenceFromMask(mask=('[#4 ]',), )
        a.Set(faces=faces1, name='Set-1')
        v1 = a.instances['Part-rigid-2'].vertices
        verts1 = v1.getSequenceFromMask(mask=('[#1 ]',), )
        a.Set(vertices=verts1, name='Set-2')
        ##??????????
        regionDef = myModel.rootAssembly.sets['Set-1']
        myModel.historyOutputRequests['H-Output-1'].setValues(variables=(
            'RF3',), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
        regionDef = myModel.rootAssembly.sets['Set-2']
        myModel.HistoryOutputRequest(name='H-Output-2',
        createStepName='Step-1', variables=('U3', ), region=regionDef,
        sectionPoints=DEFAULT, rebar=EXCLUDE)
        # Interactionģ��
        ##Tie
        s1 = a.instances['Part-concrete-1'].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#8 ]',), )
        region1 = a.Surface(side1Faces=side1Faces1, name='m_Surf-10')

        s1 = a.instances['Part-rigid-1'].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#2 ]',), )
        region2 = a.Surface(side1Faces=side1Faces1, name='s_Surf-10')
        myModel.Tie(name='Constraint-1', master=region1, slave=region2,
                              positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON,
                              thickness=ON)
        s1 = a.instances['Part-concrete-1'].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#4 ]',), )
        region1 = a.Surface(side1Faces=side1Faces1, name='m_Surf-12')
        s1 = a.instances['Part-rigid-2'].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#4 ]',), )
        region2 = a.Surface(side1Faces=side1Faces1, name='s_Surf-12')
        myModel.Tie(name='Constraint-2', master=region1, slave=region2,
                              positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON,
                              thickness=ON)


        myModel.ContactProperty('IntProp-1')
        myModel.interactionProperties['IntProp-1'].NormalBehavior(
            pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)
        myModel.interactionProperties['IntProp-1'].TangentialBehavior(
            formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
            pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=(( 0.25,),), shearStressLimit=None,
            maximumElasticSlip=FRACTION,
            fraction=0.005, elasticSlipStiffness=None)
        s1 = a.instances['Part-concrete-1'].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#1 ]',), )
        region1 = a.Surface(side1Faces=side1Faces1, name='m_Surf-13')

        s1 = a.instances['Part-tube-1-out'].faces
        side2Faces1 = s1.getSequenceFromMask(mask=('[#1 ]',), )
        region2 = a.Surface(side2Faces=side2Faces1, name='s_Surf-13')
        myModel.SurfaceToSurfaceContactStd(name='Int-1',
                                                     createStepName='Initial', master=region1, slave=region2, sliding=FINITE,
                                                     thickness=ON, interactionProperty='IntProp-1', adjustMethod=NONE,
                                                     initialClearance=OMIT, datumAxis=None, clearanceRegion=None)
        s1 = a.instances['Part-tube-2-in'].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#1 ]',), )
        region1 = a.Surface(side1Faces=side1Faces1, name='m_Surf-7')

        s1 = a.instances['Part-concrete-1'].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#2 ]',), )
        region2 = a.Surface(side1Faces=side1Faces1, name='s_Surf-7')
        myModel.SurfaceToSurfaceContactStd(name='Int-2',
                                                     createStepName='Initial', master=region1, slave=region2,
                                                     sliding=FINITE, thickness=ON, interactionProperty='IntProp-1',
                                                     adjustMethod=NONE, initialClearance=OMIT, datumAxis=None,
                                                     clearanceRegion=None)

        s1 = a.instances['Part-tube-1-out'].edges
        side1Edges1 = s1.getSequenceFromMask(mask=('[#1 ]',), )
        s2 = a.instances['Part-tube-2-in'].edges
        side1Edges2 = s2.getSequenceFromMask(mask=('[#1 ]',), )
        region1 = a.Surface(side1Edges=side1Edges1 + side1Edges2, name='m_Surf-9')
        s1 = a.instances['Part-rigid-2'].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#4 ]',), )
        region2 = a.Surface(side1Faces=side1Faces1, name='s_Surf-9')
        myModel.ShellSolidCoupling(name='Constraint-3', shellEdge=region1,
                                             solidFace=region2, positionToleranceMethod=COMPUTED)


        s1 = a.instances['Part-tube-2-in'].edges
        side1Edges1 = s1.getSequenceFromMask(mask=('[#2 ]',), )
        s2 = a.instances['Part-tube-1-out'].edges
        side1Edges2 = s2.getSequenceFromMask(mask=('[#2 ]',), )
        region1 = a.Surface(side1Edges=side1Edges1 + side1Edges2, name='m_Surf-11')

        s1 = a.instances['Part-rigid-1'].faces
        side1Faces1 = s1.getSequenceFromMask(mask=('[#2 ]',), )
        region2 = a.Surface(side1Faces=side1Faces1, name='s_Surf-11')
        myModel.ShellSolidCoupling(name='Constraint-4', shellEdge=region1,
                                             solidFace=region2, positionToleranceMethod=COMPUTED)

        # Load???
        f1 = a.instances['Part-rigid-1'].faces
        faces1 = f1.getSequenceFromMask(mask=('[#4 ]',), )
        region = regionToolset.Region(faces=faces1)
        myModel.DisplacementBC(name='BC-1', createStepName='Initial',
                               region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                               amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
        f1 = a.instances['Part-rigid-2'].faces
        faces1 = f1.getSequenceFromMask(mask=('[#2 ]',), )
        region = regionToolset.Region(faces=faces1)
        myModel.DisplacementBC(name='BC-2', createStepName='Initial',
                               region=region, u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
                               amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
        myModel.boundaryConditions['BC-2'].setValuesInStep(stepName='Step-1', u3=(-1 * displacement))

        # Mesh
        a = myModel.rootAssembly
        e1 = a.instances['Part-rigid-2'].edges
        v1 = a.instances['Part-rigid-2'].vertices
        vv2 = v1.findAt((diameter_out/ 2, 0, length), )
        vv3 = v1.findAt((diameter_out / 2, 0, length + 20), )
        DP = a.DatumPlaneByThreePoints(point2=vv2, point3=vv3,
                                       point1=a.instances['Part-rigid-2'].InterestingPoint(edge=e1[0], rule=MIDDLE))

        c1 = a.instances['Part-concrete-1'].cells
        cells1 = c1.getSequenceFromMask(mask=('[#1 ]',), )
        c2 = a.instances['Part-rigid-1'].cells
        cells2 = c2.getSequenceFromMask(mask=('[#1 ]',), )
        c3 = a.instances['Part-rigid-2'].cells
        cells3 = c3.getSequenceFromMask(mask=('[#1 ]',), )
        pickedCells = cells1 + cells2 + cells3
        d4 = a.datums
        a.PartitionCellByDatumPlane(datumPlane=d4[DP.id], cells=pickedCells)
        e1 = a.instances['Part-rigid-2'].edges
        e2 = a.instances['Part-rigid-1'].edges

        f1 = a.instances['Part-tube-1-out'].faces
        faces1 = f1.getSequenceFromMask(mask=('[#1 ]',), )
        f2 = a.instances['Part-tube-2-in'].faces
        faces2 = f2.getSequenceFromMask(mask=('[#1 ]',), )
        pickedFaces = faces1 + faces2
        d1 = a.datums
        a.PartitionFaceByDatumPlane(datumPlane=d1[DP.id], faces=pickedFaces)
        e1 = a.instances['Part-rigid-2'].edges
        e2 = a.instances['Part-rigid-1'].edges
        
        
        DP = a.DatumPlaneByThreePoints(point1=a.instances['Part-rigid-2'].InterestingPoint(edge=e1[4], rule=MIDDLE),
                                       point2=a.instances['Part-rigid-2'].InterestingPoint(edge=e1[7], rule=MIDDLE),
                                       point3=a.instances['Part-rigid-2'].InterestingPoint(edge=e1[5], rule=MIDDLE))
        c1 = a.instances['Part-concrete-1'].cells
        cells1 = c1.getSequenceFromMask(mask=('[#3 ]',), )
        c2 = a.instances['Part-rigid-1'].cells
        cells2 = c2.getSequenceFromMask(mask=('[#3 ]',), )
        c3 = a.instances['Part-rigid-2'].cells
        cells3 = c3.getSequenceFromMask(mask=('[#3 ]',), )
        pickedCells = cells1 + cells2 + cells3
        d4 = a.datums
        a.PartitionCellByDatumPlane(datumPlane=d4[DP.id], cells=pickedCells)
        ###��Y-Z�������и����
        f1 = a.instances['Part-tube-1-out'].faces
        faces1 = f1.getSequenceFromMask(mask=('[#3 ]',), )
        f2 = a.instances['Part-tube-2-in'].faces
        faces2 = f2.getSequenceFromMask(mask=('[#3 ]',), )
        pickedFaces = faces1 + faces2
        d1 = a.datums
        a.PartitionFaceByDatumPlane(datumPlane=d1[DP.id], faces=pickedFaces)  # ����DatumPlane��id���ж�λ
        ###�ṹ������




        f1 = a.instances['Part-tube-1-out'].faces
        faces1 = f1.getSequenceFromMask(mask=('[#f ]',), )
        pickedRegions = faces1
        a.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
        f2 = a.instances['Part-tube-2-in'].faces
        faces2 = f2.getSequenceFromMask(mask=('[#f ]',), )
        pickedRegions = faces2
        a.setMeshControls(regions=pickedRegions, technique=STRUCTURED)


###正式划分网格

        a = myModel.rootAssembly
        e1 = a.instances['Part-tube-1-out'].edges
        e2 = a.instances['Part-tube-2-in'].edges
        pickedEdges = e1.getSequenceFromMask(mask=('[#eea ]',), ) + \
                      e2.getSequenceFromMask(mask=('[#eea ]',), )
        a.seedEdgeByNumber(edges=pickedEdges, number=meshsize, constraint=FIXED)
        
        a = myModel.rootAssembly
        e1 = a.instances['Part-concrete-1'].edges
        pickedEdges = e1.getSequenceFromMask(mask=('[#f1f56d00 ]',), )
        a.seedEdgeByNumber(edges=pickedEdges, number=meshsize, constraint=FIXED)
        
        a = myModel.rootAssembly
        e1 = a.instances['Part-concrete-1'].edges
        pickedEdges = e1.getSequenceFromMask(mask=('[#20a02aa ]',), )
        a.seedEdgeByNumber(edges=pickedEdges, number=5, constraint=FIXED)
        
        a = myModel.rootAssembly
        partInstances = (a.instances['Part-concrete-1'],
                         a.instances['Part-tube-1-out'], a.instances['Part-tube-2-in'],)
        a.seedPartInstance(regions=partInstances, size=25.0, deviationFactor=0.1,
                            minSizeFactor=0.1)
        a = myModel.rootAssembly
        partInstances = (a.instances['Part-concrete-1'],
                         a.instances['Part-tube-1-out'], a.instances['Part-tube-2-in'],)
        a.generateMesh(regions=partInstances)
        a = myModel.rootAssembly
        partInstances = (a.instances['Part-rigid-1'], a.instances['Part-rigid-2'],)
        a.seedPartInstance(regions=partInstances, size=meshsize*4, deviationFactor=0.1,
                           minSizeFactor=0.1)
        a = myModel.rootAssembly
        partInstances = (a.instances['Part-rigid-1'], a.instances['Part-rigid-2'],)
        a.generateMesh(regions=partInstances)


        viewportshow.assemblyDisplay.setValues(renderStyle=SHADED, mesh=ON)
        viewportshow.setValues(displayedObject=a)

        mdb.Job(name=jobname, model=Modelname, description='', type=ANALYSIS,
                    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=50,
                    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
                    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
                    scratch='', multiprocessingMode=DEFAULT, numCpus=1)
        # ??????????
        session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
        reply = getWarningReply(message='Would you like to submit the job?', buttons=(YES, NO))
        if reply == YES:
            mdb.jobs[jobname].submit(consistencyChecking=OFF)
