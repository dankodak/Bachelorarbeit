Für das Greedy-Verfahren:
setting = 'weighted'|'standard'
	gewichtetes oder Standardverfahren
n = Anzahl der Iterationen
test_grid = 0|1
	Validationspunkte zufällig oder auf Gitter
test = Anzahl der Validationspunkte
symmetric = 0|1
	Nicht-Symmetrisches oder Symmetrisches Verfahren
kernel = 'gauss'|'wendland'
pde = 'square'|'circle'|'disc'
error = 'abs'|'res'
	Fehlerausgabe absolut oder im Residuum

-------------------------------------------------------------

Zufällige Kollokationspunkte oder auf Gitter:
setting = 'weighted'|'standard'
	gewichtetes oder Standardverfahren
grid = 0|1
	Kollokationspunkte zufällig oder auf Gitter
test_grid = 0|1
	Validationspunkte zufällig oder auf Gitter
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

Bei über einem Gitter gewählten Punkten (Validations- und Kollokationspunkte) wird versucht ein nxn Gitter über das Gebiet zu legen. Man hat also < n^2 Kollokationspunkte.

In den solvePDE Dateien kann eine Konditionsberechnung und deren Plot auskommentiert werden.

Zur Definition einer PDE muss sie bei allFunctions angelegt werden. Dort muss im switch pde ein neuer case angelegt werden mit:
	-ws = Gewichtsfunktion
	-fs = Rechte Seite der PDE
	-realSols = 0|echte Lösung. Wird 0 eingegeben, wird ein Fehlerschätzer im Residuum benutzt, wird eine Lösung eingegeben, ein absoluter.
	-realSolPlots = Echte Lösung für Fehlerberechnung auf Testpunkten.
Zur Definition eines neuen Kerns muss im switch kernel ein neuer case angelegt werden. Dabei ist rbfs der Kern.
lap_rbf und lap2_rbf müssen dann jeweils der auf rbfs angewandte Differentialoperator sein.

Für die Standardverfahren müssen für jedes Gebiet noch die Randpunkte in boundary_points.m von Hand programmiert werden.


--------------------------------------------------------------

Die 3D und 4D Versionen sind von der Bedienung her gleich, aber eingeschränkt. Die beiden Programme wurden nur dafür benutzt, die in der Arbeit gezeigten PDEs zu lösen.
Sie sind von daher nur auf eben diese PDEs angepasst und können ohne weiteres keine anderen lösen (Hauptsächlich, weil die Kollokationspunkte nicht sauber umgeschrieben wurden).
Ebenso beherrschen sie deswegen auch nur die gewichtete Kollokation.