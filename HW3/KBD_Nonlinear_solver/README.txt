GENERAL FINITE ELEMENT LIBRARY (+fea)
====================================
The package includes Q4 and Q8 plane-strain elements, shared kinematics,
Total Lagrangian element response, material models, global assembly, solver
classes, reference edge tractions, mesh import and contour plotting.

Files:
+fea/+elements/Quad4.m, Quad8.m         Element-specific interpolation
+fea/+elements/ElementBase.m            Shared public element API
+fea/+elements/ElementCalculations.m    All shared TL/strain/integration math
+fea/+materials/                       StVK and alternative log material
+fea/+solvers/                          Linear, Newton, internal penalty helper
+fea/+loads/SurfaceTraction.m            Fixed reference-edge traction
+fea/+assembly/Assembler.m              Q4/Q8 global assembly
+fea/+mesh/AbaqusINP.m                  Existing mesh reader (Q4-focused)
+fea/+results/FieldPlotter.m            Existing field contour plotter

Use from the directory containing +fea:
mat=fea.materials.LinearElastic(70000,0.3);
q4=fea.elements.Quad4([0 0;1 0;1 1;0 1],mat,1);
[fint,KL,KNL]=q4.Element_Response_cal(zeros(8,1));

For Q8 use nodes [corner1..corner4, midpoint12,23,34,41]
and 3x3 integration. Existing Abaqus mesh reader has not been
upgraded to import Q8 mesh; Q8 analysis accepts supplied coordinates
and 8-node connectivity via Assembler.Assemble_TL_cal.

Migration: replace fea_hw2. with fea. in existing driver scripts.
Old calls to TotalLagrangian.Element_Response_cal(element,mat,q) become
 element.Element_Response_cal(q) or
 fea.elements.ElementCalculations.Element_Response_cal(element,mat,q).
Old FiniteStrain methods now live in fea.elements.ElementCalculations.
Do NOT remove your old package until all original tests pass.

Notes: MATLAB execution was unavailable in this environment; run
examples/Verify_Quad4_Quad8.m then regression-test the existing cases.
