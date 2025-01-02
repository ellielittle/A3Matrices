load "ft_main.m";

Q := Integers();

GraphAuto := function(x)
    i := 1;
    while i le #x do
        if x[i] eq 3+1 
            then x[i] := 1;
        else 
            x[i] := x[i]+1;
        end if;
        i:=i+1;
        end while;
    return x;
end function;

GraphAutoN := function(x, j)
    i := 1;
    while i le j do 
        x := GraphAuto(x);
        i := i+1;
    end while;
    return x;
end function;

matrixrowdata := function(u,t,datum,Jdatum,var)
    useq := Eltseq(Jdatum[3][u]);
    pos_root := [Transpose(Matrix(r)) : r in PositiveRoots(datum[3])];
    if t eq [] then 
        return [<1,Jdatum[3][u]>];
    else
        a := datum[4][t[1]];
            if CheckInRoots(useq, a, Jdatum[2], datum) then
                return [<-1/var[1], Jdatum[3][u]>];
            elif CheckInRoots(useq, a, pos_root, datum) then
                return [<(var[1]-(1/var[1])),datum[1] ! useq>,<1,(datum[1] ! useq)*(datum[1] ! t)>];
            else
                return [<1,(datum[1] ! useq)*(datum[1] ! t)>];
            end if;
    end if;
    end function;

WJreduce := function(w, W, Jdatum)
    J := Jdatum[1];
    Jelts := [Eltseq(x) : x in Jdatum[3]];
    Jelts := [W ! x : x in Jelts];
    if W ! w in Jelts then 
        return w;
    else
        neww := w;
        while (W ! neww) in Jelts eq false do 
            for i in J do
                if Length(W ! ([i] cat neww)) le Length(W ! neww) then
                    neww := Eltseq(W ! ([i] cat neww));
                end if;
            end for;
        end while;
        return neww;
    end if;
    end function;

WJcomp := function(t,datum,Jdatum);
    J := Jdatum[1];
    W := datum[1];
    theta := transDecomp(t[1],datum)[1];
    M := Matrix([[2/1,-1/1,0/1],[-1/1,2/1,-1/1],[0/1,-1/1,2/1]]);
    wt := Eltseq(M*Transpose(Matrix([[x/1 : x in transDecomp(t[1],datum)[2]]])));
    finite := GraphAutoN(Eltseq(theta), 4-t[2]);

    if t[2] eq 1 then 
        new := <[1,2,3],[1/1,0,0]>;
    elif t[2] eq 2 then 
        new := <[2,1,3,2],[0,1/1,0]>;
    elif t[2] eq 3 then 
        new := <[3,2,1],[0,0,1/1]>;
    elif t[2] eq 0 then 
        new := <[],[0,0,0]>;
    else 
        return "error";
    end if;

    neweight := Eltseq(M*Transpose(Matrix([[x/1 : x in transDecomp(Eltseq(W ! (new[1] cat finite)),datum)[2]]])));

    wt := [wt[i] + new[2][i] + neweight[i]/1 : i in [1..3]];
    finite := transDecomp(Eltseq(W ! (new[1] cat finite)),datum)[1];
    Ematrix := Matrix([[1/1,1/1,1/1],[0/1,1/1,1/1],[0/1,0/1,1/1],[0/1,0/1,0/1]]);
    wte := Eltseq(Ematrix * Transpose(Matrix([wt])));
    return WJreduce(Eltseq(finite),W,Jdatum),wte;
    end function;

zpart := function(wt,var,J)
    z := [var[x] : x in [2..#var]];
    zpart := 1;
    if J eq {1} or J eq {1,3} then 
        for x in [1..#wt] do 
            if x eq 1 or x eq 2 then 
                if wt[x] lt 0 then 
                    zpart := zpart*(1/z[1])^(-wt[x]);
                elif wt[x] gt 0 then 
                    zpart := zpart*z[1]^(wt[x]);
                else
                    zpart := zpart;
                end if;
            elif x eq 3 or x eq 4 then 
                if wt[x] lt 0 then 
                    zpart := zpart*(1/z[2])^(-wt[x]);
                elif wt[x] gt 0 then 
                    zpart := zpart*z[2]^(wt[x]);
                else
                    zpart := zpart;
                end if;
            else
                return "error";
            end if;
        end for;
    elif J eq {1,2} then 
        for x in [1..#wt] do 
            if wt[x] lt 0 then 
                zpart := zpart*(1/z[1])^(-wt[x]);
            elif wt[x] gt 0 then 
                zpart := zpart*z[1]^(wt[x]);
            else
                zpart := zpart;
            end if;
        end for;
    elif J eq {} then 
        for x in [1..#wt] do 
            if x eq 1 then 
                if wt[x] lt 0 then 
                    zpart := zpart*(1/z[1])^(-wt[x]);
                elif wt[x] gt 0 then 
                    zpart := zpart*z[1]^(wt[x]);
                else
                    zpart := zpart;
                end if;
            elif x eq 2  then 
                if wt[x] lt 0 then 
                    zpart := zpart*(1/z[2])^(-wt[x]);
                elif wt[x] gt 0 then 
                    zpart := zpart*z[2]^(wt[x]);
                else
                    zpart := zpart;
                end if;
            elif x eq 3 or x eq 4 then 
                if wt[x] lt 0 then 
                    zpart := zpart*(1/z[3])^(-wt[x]);
                elif wt[x] gt 0 then 
                    zpart := zpart*z[3]^(wt[x]);
                else
                    zpart := zpart;
                end if;
            else
                return "error";
            end if;
        end for;
    elif J eq {1,2,3} then 
        for x in [1..#wt] do 
            if wt[x] lt 0 then 
                zpart := zpart*(1/z[1])^(-wt[x]);
            elif wt[x] gt 0 then 
                zpart := zpart*z[1]^(wt[x]);
            else
                zpart := zpart;
            end if;
        end for;
    else 
        return "error";
    end if;
    return zpart;
    end function;

matrixrownewdata := function(u,t,datum,Jdatum,var)
    data := matrixrowdata(u,t[1],datum,Jdatum,var);
    newdata := [];
    for x in data do 
        v,wt := WJcomp(<x[2],t[2]>, datum, Jdatum);
        newdata := newdata cat [<x[1]*zpart(wt,var,Jdatum[1]),v>];
    end for;
    return newdata;
    end function;

// MIGHT NEED TO CHANGE TO INCLUDE PATHS ENDING AT SAME SPOT
rowdata := function(u,t,datum,Jdatum,var)
    data := matrixrownewdata(u,t,datum,Jdatum,var);
    datatheta := [x[2]: x in data];
    Jelts := [Eltseq(x) : x in Jdatum[3]];
    newdataa := [];
    for u in Jelts do 
        if u in datatheta then 
            newdataa := newdataa cat [data[Index(datatheta, u)]];
        else 
            newdataa := newdataa cat [<0,u>];
        end if;
    end for;
    if not #Jdatum[3] eq #newdataa then
        return "error";
    end if;
    return newdataa;
    end function;

matrixrow := function(u,t,datum,Jdatum,var)
    W := datum[1];
    data := rowdata(u,t,datum,Jdatum,var);
    data := [<x[1],W ! x[2]> : x in data];
    row := [x[1] : x in data];
    return row;
    end function;

T := function(t,datum,Jdatum,var)
    R := [];
    for u in [1..#Jdatum[3]] do
        R := R cat [matrixrow(u,t,datum,Jdatum,var)];
    end for;
    return Matrix(R);
end function;

printmatrixandleading := function(x,J)
    type := "A3";
    datum := CoxeterSetup(type);
    w := <Eltseq(datum[1] ! Eltseq(x[1])),x[2]>;
    Jdatum := ParabolicSetup(J, datum, 3);
    j := #Jdatum[3];
    if Jdatum[1] eq {1} then 
        k := 3;
        LL := LaurentSeriesRing(Q);
        AssignNames(~LL, ["z1"]);
        z1 := LL.1;
        LLL := LaurentSeriesRing(LL);
        AssignNames(~LLL, ["z2"]);
        z2 := LLL.1;
        L := LaurentSeriesRing(LLL);
        AssignNames(~L, ["q"]);
        q := L.1;
        var := [q,z1,z2];
        rows := Rows(IdentityMatrix(Q,j));
        M := Matrix([rows[5]] cat [rows[i] : i in [1..4]] cat[rows[i] : i in [6..j]]);
    elif Jdatum[1] eq {1,2} then 
        k := 1;
        LL := LaurentSeriesRing(Q);
        AssignNames(~LL, ["z1"]);
        z1 := LL.1;
        L := LaurentSeriesRing(LL);
        AssignNames(~L, ["q"]);
        q := L.1;
        var := [q,z1];
        rows := Rows(IdentityMatrix(Q,j));
        M := Matrix([rows[3]] cat [rows[1]] cat [rows[2]] cat[rows[4]]);
    elif Jdatum[1] eq {1,3} then 
        k := 2;
        LLL :=LaurentSeriesRing(Q);
        AssignNames(~LLL, ["z1"]);
        z1 := LLL.1;
        LL := LaurentSeriesRing(LLL);
        AssignNames(~LL, ["z2"]);
        z2 := LL.1;
        L := LaurentSeriesRing(LL);
        AssignNames(~L, ["q"]);
        q := L.1;
        var := [q,z1,z2];
        // Changing the ordering of {^J}W
        M := Matrix([[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]);
    elif Jdatum[1] eq {} then 
        k := 6;
        LLL := LaurentSeriesRing(Q);
        AssignNames(~LLL, ["z1"]);
        z1 := LLL.1;
        LL := LaurentSeriesRing(Q);
        AssignNames(~LL, ["z2"]);
        z2 := LL.1;
        LLLL := LaurentSeriesRing(Q);
        AssignNames(~LLLL, ["z3"]);
        z3 := LLLL.1;
        L := LaurentSeriesRing(LLLL);
        AssignNames(~L, ["q"]);
        q := L.1;
        var := [q,z1,z2,z3];
        M := IdentityMatrix(Q,j);
    elif Jdatum[1] eq {1,2,3} then 
        k := 1;
        LL := LaurentSeriesRing(Q);
        AssignNames(~LL, ["z1"]);
        z1 := LL.1;
        L := LaurentSeriesRing(LL);
        AssignNames(~L, ["q"]);
        q := L.1;
        var := [q,z1];
        M := IdentityMatrix(Q,j);
    else 
        return "error";
    end if;
    mat := [T(<[i],0>,datum,Jdatum,var) : i in [1..4]];
    X := T(<[],1>,datum,Jdatum,var);
    mat := [M*y*M: y in mat];
    X := M*X*M;
    w;
    if #w[1] eq 0 then 
        w := IdentityMatrix(L,j)*X^w[2];
    else 
        w := &*[mat[i] : i in w[1]]*X^w[2];
    end if;
    var := [x : x in var | not x eq 0];
    res := (&*[x : x in var])*Matrix([[0 : x in [1..j]] : y in [1..j]]);
    for x in [1..j] do
        for y in [1..j] do 
            res[x,y] := Coefficient(w[x,y],k);
        end for;
    end for;

    return w,res;
end function;
