F�r das Greedy-Verfahren:
setting = 'weighted'|'standard'
	gewichtetes oder Standardverfahren
n = Anzahl der Iterationen
test_grid = 0|1
	Validationspunkte zuf�llig oder auf Gitter
test = Anzahl der Validationspunkte
symmetric = 0|1
	Nicht-Symmetrisches oder Symmetrisches Verfahren
kernel = 'gauss'|'wendland'
pde = 'square'|'circle'|'disc'
error = 'abs'|'res'
	Fehlerausgabe absolut oder im Residuum

-------------------------------------------------------------

Zuf�llige Kollokationspunkte oder auf Gitter:
setting = 'weighted'|'standard'
	gewichtetes oder Standardverfahren
grid = 0|1
	Kollokationspunkte zuf�llig oder auf Gitter
test_grid = 0|1
	Validationspunkte zuf�llig oder auf Gitter
m = start:step:end
	Anzahl der Kollokationspunkte
test = Anzahl der Validationspunkte
symmetric = 0|1
	Nicht-Symmetrisches oder Symmetrisches Verfahren
kernel = 'gauss'|'wendland'
pde = 'square'|'circle'|'disc'
error = 'abs'|'res'
	Fehlerausgabe absolut oder im Residuum

--------------------------------------------------------------

Bei �ber einem Gitter gew�hlten Punkten (Validations- und Kollokationspunkte) wird versucht ein nxn Gitter �ber das Gebiet zu legen. Man hat also < n^2 Kollokationspunkte.

In den solvePDE Dateien kann eine Konditionsberechnung und deren Plot auskommentiert werden.

Zur Definition einer PDE muss sie bei allFunctions angelegt werden. Dort muss im switch pde ein neuer case angelegt werden mit:
	-ws = Gewichtsfunktion
	-fs = Rechte Seite der PDE
	-realSols = 0|echte L�sung. Wird 0 eingegeben, wird ein Fehlersch�tzer im Residuum benutzt, wird eine L�sung eingegeben, ein absoluter.
	-realSolPlots = Echte L�sung f�r Fehlerberechnung auf Testpunkten.
Zur Definition eines neuen Kerns muss im switch kernel ein neuer case angelegt werden. Dabei ist rbfs der Kern.
lap_rbf und lap2_rbf m�ssen dann jeweils der auf rbfs angewandte Differentialoperator sein.

F�r die Standardverfahren m�ssen f�r jedes Gebiet noch die Randpunkte in boundary_points.m von Hand programmiert werden.


--------------------------------------------------------------

Die 3D und 4D Versionen sind von der Bedienung her gleich, aber eingeschr�nkt. Die beiden Programme wurden nur daf�r benutzt, die in der Arbeit gezeigten PDEs zu l�sen.
Sie sind von daher nur auf eben diese PDEs angepasst und k�nnen ohne weiteres keine anderen l�sen (Haupts�chlich, weil die Kollokationspunkte nicht sauber umgeschrieben wurden).
Ebenso beherrschen sie deswegen auch nur die gewichtete Kollokation.