classdef Assembler
    % Assembly for 2D Q4 and Q8 continuum elements (2 DOFs per node).
    methods (Static)
        function [fint,KL,KNL,KT]=Assemble_TL_cal(nodes,conn,material,thickness,qGlobal,elementType)
            if nargin<6 || isempty(elementType)
                if size(conn,2)==4,elementType='Quad4';
                elseif size(conn,2)==8,elementType='Quad8';
                else,error('Unsupported element connectivity');end
            end
            nNodes=size(nodes,1);nDOF=2*nNodes;
            assert(numel(qGlobal)==nDOF,'Global DOF mismatch');
            assert(all(conn(:)>=1 & conn(:)<=nNodes),'Invalid connectivity');
            fint=zeros(nDOF,1);KL=zeros(nDOF);KNL=zeros(nDOF);
            for e=1:size(conn,1)
                ids=conn(e,:);assert(numel(unique(ids))==numel(ids),'Repeated element nodes');
                dofs=reshape([2*ids-1;2*ids],[],1);
                switch lower(elementType)
                    case {'quad4','cpe4'},el=fea.elements.Quad4(nodes(ids,:),material,thickness);
                    case {'quad8','cpe8'},el=fea.elements.Quad8(nodes(ids,:),material,thickness);
                    otherwise,error('Unknown element type');
                end
                [fe,ke,kg]=el.Element_Response_cal(qGlobal(dofs));
                fint(dofs)=fint(dofs)+fe;
                KL(dofs,dofs)=KL(dofs,dofs)+ke;
                KNL(dofs,dofs)=KNL(dofs,dofs)+kg;
            end
            KT=KL+KNL;
        end
        function [fint,KL,KNL,KT]=Assemble_TL_Q4_cal(nodes,conn,material,thickness,qGlobal,~)
            [fint,KL,KNL,KT]=fea.assembly.Assembler.Assemble_TL_cal( ...
                nodes,conn,material,thickness,qGlobal,'Quad4');
        end
    end
end
