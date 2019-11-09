import sys
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5 import uic
from BlackScholesPDECalculatorCode import EuropeanOptions

qtCreatorFile = "Black Scholes PDE Calculator.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class Main(QMainWindow, Ui_MainWindow):          
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.ExplicitFDM.setChecked(True)
        self.Calculate.clicked.connect(self.OptionPricing)
        
    def OptionPricing(self):
        S = float(self.SharePriceInput.text())
        K = float(self.StrikePriceInput.text())
        r = float(self.RiskFreeInterestRateInput.text())
        q = float(self.SharesDividendYieldInput.text())
        T = float(self.OptionsTimetoMaturityInput.text())
        sigma = float(self.VolatilityofSharePriceInput.text())
        M = int(self.StepSizeforSharePriceMInput.text())
        N = int(self.StepSizeforTimeNInput.text())           
        V = EuropeanOptions(S, K, r, q, T, sigma)
        
        if self.ExplicitFDM.isChecked():
            Value = V.ExplicitFDM(M, N)
        elif self.ImplicitFDM.isChecked():
            Value = V.ImplicitFDM(M, N)
        else:
            Value = V.CrankNicolsonM(M, N)
            
        self.CallOutput.setText(str(round(Value['call'],4)))
        self.PutOutput.setText(str(round(Value['put'],4)))

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())