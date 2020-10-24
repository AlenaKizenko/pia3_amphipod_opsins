// main directory
mdir = "D:\QE_Pro";

// current subdirectories in the main one
sdirs = ["2019.11.22 LED in box\blue"
         "2019.11.22 LED in box\green"
         "2019.11.22 LED in box\yellow"
         "2019.11.22 LED in box\red"
         ];

// names of files in corresponding subdirectories
fin1 = [];
fin2 = [];
fin3 = [];
fin4 = [];
fin5 = [];
fin6 = [];
fin7 = [];
fin8 = [];
fin9 = [];
fin10 = [];
fin11 = [];
fin12 = [];
fin13 = [];
fin14 = [];
fin15 = [];
fin16 = [];
fin17 = [];
fin18 = [];
fin19 = [];
if ( (sum(size(sdirs)) > 1) & (sum(length(fin1)) == 0) ) then fin1 = (basename(listfiles(fullfile(mdir, sdirs(1), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 2) & (sum(length(fin2)) == 0) ) then fin2 = (basename(listfiles(fullfile(mdir, sdirs(2), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 3) & (sum(length(fin3)) == 0) ) then fin3 = (basename(listfiles(fullfile(mdir, sdirs(3), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 4) & (sum(length(fin4)) == 0) ) then fin4 = (basename(listfiles(fullfile(mdir, sdirs(4), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 5) & (sum(length(fin5)) == 0) ) then fin5 = (basename(listfiles(fullfile(mdir, sdirs(5), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 6) & (sum(length(fin6)) == 0) ) then fin6 = (basename(listfiles(fullfile(mdir, sdirs(6), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 7) & (sum(length(fin7)) == 0) ) then fin7 = (basename(listfiles(fullfile(mdir, sdirs(7), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 8) & (sum(length(fin8)) == 0) ) then fin8 = (basename(listfiles(fullfile(mdir, sdirs(8), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 9) & (sum(length(fin9)) == 0) ) then fin9 = (basename(listfiles(fullfile(mdir, sdirs(9), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 10) & (sum(length(fin10)) == 0) ) then fin10 = (basename(listfiles(fullfile(mdir, sdirs(10), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 11) & (sum(length(fin11)) == 0) ) then fin11 = (basename(listfiles(fullfile(mdir, sdirs(11), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 12) & (sum(length(fin12)) == 0) ) then fin12 = (basename(listfiles(fullfile(mdir, sdirs(12), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 13) & (sum(length(fin13)) == 0) ) then fin13 = (basename(listfiles(fullfile(mdir, sdirs(13), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 14) & (sum(length(fin14)) == 0) ) then fin14 = (basename(listfiles(fullfile(mdir, sdirs(14), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 15) & (sum(length(fin15)) == 0) ) then fin15 = (basename(listfiles(fullfile(mdir, sdirs(15), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 16) & (sum(length(fin16)) == 0) ) then fin16 = (basename(listfiles(fullfile(mdir, sdirs(16), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 17) & (sum(length(fin17)) == 0) ) then fin17 = (basename(listfiles(fullfile(mdir, sdirs(17), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 18) & (sum(length(fin18)) == 0) ) then fin18 = (basename(listfiles(fullfile(mdir, sdirs(18), "*.txt"))))' + ".txt"; end
if ( (sum(size(sdirs)) > 19) & (sum(length(fin19)) == 0) ) then fin19 = (basename(listfiles(fullfile(mdir, sdirs(19), "*.txt"))))' + ".txt"; end
fl = tlist(fin1, fin2, fin3, fin4, fin5, fin6, fin7, fin8, fin9, fin10, fin11, fin12, fin13, fin14, fin15, fin16, fin17, fin18, fin19);

// transform the list into matrix
n = 0
for i = 1:length(fl)
    finl = size(fl(i));
    n = max(n, finl(2));
    if (sum(size(fl(i)) > 0)) > 0 then
        m = i;
    end
end
fm = zeros(m, n)
fm = string(fm)
fm(:,:) = ""
l = []
for i = 1:m
    finl = size(fl(i));
    l(i) = finl(2);
    fm(i, 1:l(i)) = fl(i);
end


// load spectra data
k = 0
X = zeros(2000, sum(l))
Z = zeros(2000, sum(l))
X(:,:) = %nan
Z(:,:) = %nan
lab = []
for i = 1:m
    for j = 1:l(i)
        k = k + 1;
        M = csvRead(fullfile(mdir,sdirs(i),fm(i,j)), "\t", ",", [], [], "/:|_|>/");
        n = size(M)
        n = n(1)
        X(1:n,k) = M(:,1);    // wavelength
        Z(1:n,k) = M(:,2);    // intensity counts
        lab(k) = strcat([string(k), "_", sdirs(i), "/", fm(i,j)])    // spectra labels
    end
end

//subtract a spectrum from a spectrum
// x - extract from, y - to extract, mult - multiply y on
// would be good to add multiple subtraction...
function sds = subs(x, y, mult)
    if length(y) < length(x) then
        y = repmat(y, 1, length(x))
    end
    if length(mult) < length(x) then
        mult = repmat(mult, 1, length(x))
    end
    if sum(y == 0) > 0 then
        nosub = (y == 0)
        y(nosub) = x(nosub)
        mult(nosub) = 0
    end
    Xz = X(:,x)
    Xy = X(:,y)
    Zc = Z(:,x)
    Yc = Z(:,y)
    sds = ones(2000, length(x))
    for i = 1:length(x)
        if sum(~( (Xz(:,i) == Xy(:,i)) | (isnan(Xz(:,i))) )) == 0 then
            sds(:,i) = Zc(:,i) - Yc(:,i)*mult(i)
        else        // this part is not checked for errors
            for j = 1:2000
                sds(j,i) = Zc(j,i) - mean(Yc( (Xy(:,i) >= Xz(j,i) - 1) & (Xy(:,i) <= Xz(j,i) + 1), i )) * mult(i)
            end
        end
    end
endfunction


// apply simple moving average to Z
// x - numbers of columns in X, mw - size of moving window
function zf = smmz(x, mw, y, mult)
    if sum(y > 0) == 0 then
        Zc = Z(:,x)
    else
        Zc = subs(x, y, mult)
    end
    Zb = Zc
    Xc = X(:,x)
    if length(mw) == 1 then
        mw = repmat(mw, 1, length(x))
    end
    for i = 1:length(x)
        for j = 1:2000
            Zc(j,i) = mean(Zb( (Xc(:,i) >= Xc(j,i) - mw(i)) & (Xc(:,i) <= Xc(j,i) + mw(i)), i ))
        end
    end
    zf = Zc;
endfunction


// calculate ratio between 2 spectra at a certain range/point
// x - two columns of X, wr - range(s) of wavelengthes
// y and mult — pass arguments to subs()
// 1 range (2 numbers) will be applied to both spectra, 2 ranges (4 numbers) - respectively
function ratio = rt(x, wr, y, mult) 
    if sum(y > 0) == 0 then
        Xc = X(:,x)
        Zc = Z(:,x)
    else
        Xc = X(:,x)
        Zc = subs(x, y, mult)
    end
    if length(wr) == 2 then
        m1 = mean( Zc( (Xc(:,1) >= wr(1)) & (Xc(:,1) <= wr(2)), 1) )
        m2 = mean( Zc( (Xc(:,2) >= wr(1)) & (Xc(:,2) <= wr(2)), 2) )
    else
        m1 = mean( Zc( (Xc(:,1) >= wr(1)) & (Xc(:,1) <= wr(2)), 1) )
        m2 = mean( Zc( (Xc(:,2) >= wr(3)) & (Xc(:,2) <= wr(4)), 2) )
    end
    ratio = m1/m2
endfunction


// function to draw
// x - columns of X, mw - moving window, mark - mark to draw,
// mult - multiply x on
function psm(x, mw, mark, mult, y, ymult)
    clf()
    if length(mult) == 1 then
        mult = repmat(mult, 1, length(x))
    end
    mult = repmat(mult, 2000, 1)
    if mw == 0 then
        Xm = X(:,x)
        Zm = subs(x, y, ymult)
        Zm = Zm.*mult
        Zf = Zm
    else
        Xm = X(:,x)
        Zm = subs(x, y, ymult)
        Zm = Zm.*mult
        Zf = smmz(x, mw, y, ymult);
        Zf = Zf.*mult
    end
    colvec = [2,3,5,4,6,7,1] // WTF with color schemes in this program?
    colvec = repmat(colvec,1,100)
    plot(Xm, Zf);
    if mark <> "" then 
        plot(Xm, Zm, mark)
    end
    plw=gca();
    plw.margins=[0.08,0.02,0.05,0.04];
    legends((lab(x))', colvec(1:length(x)), opt=6)
    xgrid()
endfunction


// fit observations to known spectra
// with regression analysis and check quality of the regression
// X must be identical!
// lr — borders of wavelength range to use
// x — observations to decompose to "known" spectra
// w — "known" spectra to use in decomposition
// const — allow constants? draw — plot residual values?
function b = reg(lr, x, w, const, draw)
    xlr = (X(:,x(1)) > lr(1)) & (X(:,x(1)) < lr(2))
    if const then
        W = [ Z(xlr, w), ones( length(X(xlr, x(1))), 1 ) ]
    else
        W = [ Z(xlr, w) ]
    end
    G = Z(xlr, x)
    sw = size(W)
    b = ones(sw(2), length(x))
    for i = 1:length(x)
        b(:,i) = regres(G(:,i), W)
    end
    
    if draw then
        if const then
            Zs = [Z(:,w), ones(2000,1)]
        else
            Zs = Z(:,w)
        end
        Zf = Z(:,x) - Zs*b
        clf()
        plot(X(:,x), Zf);
        colvec = [2,3,5,4,6,7,1] // WTF with color schemes in this program?
        colvec = repmat(colvec,1,100)
        plw=gca();
        plw.margins=[0.08,0.02,0.05,0.04];
        legends((lab(x))', colvec(1:length(x)), opt=6)
        xgrid()
    end
endfunction

B = zeros(3, sum(l))
B(:, 4:sum(l)) = reg([599.5, 650], x = [4:sum(l)], w = [1:3], const = %F, draw = %F)


// plot multiple spectra with subtraction and return ratio within the same spectrum
// x — spectra, y — to substract, wr — wavelength ranges, draw — plot?
// nb — number of raw in B
// subtract just ONE spectrum!
function mrt = pmrt(x, y, wr, nb, draw)
    if draw then
        psm(x, 0, ".", 1, y, ymult = B(nb,x))
    end
    mrt = zeros(1, length(x))
    if (length(y) < length(x)) then
        y = repmat(y, length(x), 1)
    end
    for i = 1:length(x)
        mrt(i) = rt( [x(i),x(i)] , wr, y(i), mult = B(nb,x(i)))
    end
endfunction


// psm with automatic scaling of all spectra to the maximal one at certain range
// x — spectra, y — to substract, wr — wavelength range to calculate rt()
// mw, makr — to pass to psm()
// avf — write average spectra to the specified file
function psmsc(x, y, wr, mw, mark, avw)
    if length(y) < length(x) then
        y = repmat(y, 1, length(x))
    end
    mx = x(1)
    for i = 2:length(x)
        mxc = mean(Z((X(:,mx)   >= wr(1)) & (X(:,mx)   <= wr(2)), mx))
        mxp = mean(Z((X(:,x(i)) >= wr(1)) & (X(:,x(i)) <= wr(2)), x(i)))
        if (mxp > mxc) then
            mx = x(i)
        end
    end
    mu = ones(1, length(x))
    for i = 1:length(x)
        mu(i) = rt([mx, x(i)], wr, y(i), 1)
    end
    psm(x=x, mw=mw, mark=mark, mult=mu, y=y, ymult=1)
    if (length(avw) > 0) then
        avXZ = repmat(X(:,1), 1, 2)
        if y == 0 then
            for i = 1:2000
                avXZ(i, 2) = mean( (Z(i,x)).*mu )
            end
        else
            for i = 1:2000
                avXZ(i, 2) = mean( (Z(i,x) - Z(i,y)).*mu )
            end
        end
        csvWrite(avXZ, fullfile(mdir, avw), ascii(9), ",", 6)
    end
endfunction


B = ones(1, sum(l))
//psm(x = 1:sum(l), mw = 0, mark = ".", mult = 1, y = 0, ymult = 1)


// in the box, with scatterer
// take measurements on the left or right side, they are similar; in the center it's higher
// background is totally absent

pmrt(x = [1 2 4 3], y = 0, wr = [604, 606, 639, 641], 1, draw = %F)


cal = csvRead(fullfile(mdir, "Photon calibration/matrix grid and scatterer.csv"))
//Z(1:1044,21) = Z(1:1044,21)./cal
Z(1:1044,1:sum(l)) = Z(1:1044,1:sum(l))./repmat(cal, 1, sum(l))

pmrt(x = [1 2 4 3], y = 0, wr = [604, 606, 639, 641], 1, draw = %T)


// from 380 to 760 nm
bl = sum(Z(43:535, 1))
gr = sum(Z(43:535, 2))
ye = sum(Z(43:535, 3))
re = sum(Z(43:535, 4))

[bl gr ye re]/10000

csvWrite(Z(1:1044,1:sum(l)), fullfile(mdir, "2019.11.22 LED in box\calibrated.csv"), ascii(9), ",")

