function[Cp,Bp,Cm,Bm]= apply_create_anihilate_onvector(C,B,i)
    Bp = []; Cp=0;
    Bm = []; Cm=0;
    Nst = size(B,2);
    irp = 0; irm = 0;
    for ir = 1:size(B,1)
        b = B(ir,:);
        % transform to 1's and 0's representation
        V1 = b; V1(V1==-1)=0; V1(V1==2)=1; V2 = b; V2(V2==1)=0;
        V2(V2==2)=1; V2=abs(V2); V = [V1,V2];

        [Vp,signop] = apply_creation_operator(i,V);
        if (signop ~=0)
            irp = irp + 1;
            V1 = Vp(1:Nst);  V2 = 0.5*Vp(Nst+1:end); bp = V1+V2; bp(bp==1.5)=2; bp(bp==0.5)=-1;
            Bp = [Bp;bp]; Cp(irp) = C(ir)*signop;
        end

        [Vm,signom] = apply_anihilation_operator(i,V);
        if (signom ~=0)
            irm = irm + 1;
            V1 = Vm(1:Nst);  V2 = 0.5*Vm(Nst+1:end); bm = V1+V2; bm(bm==1.5)=2; bm(bm==0.5)=-1;
            Bm = [Bm;bm]; Cm(irm) = C(ir)*signom;
        end
   
    end % ir
    Cp = Cp';
    Cm = Cm';
end % end function