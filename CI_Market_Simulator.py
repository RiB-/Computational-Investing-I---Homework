# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 10:52:47 2015

@author: Rfoti
"""
#******************************************************************************
# Importing packages
#******************************************************************************
import numpy as np  #library for matrix-array analysis
import pandas as pd  #library for advanced data analysis
import pandas.io.parsers as pd_par #import parser
import matplotlib.pyplot as plt #import plotting library
import copy #import copying module
import datetime as dt #import module for date/time
#-------------------
# Importing QSTK modules
#-------------------
import QSTK.qstkutil.qsdateutil as du
import QSTK.qstkutil.DataAccess as da
import QSTK.qstkstudy.EventProfiler as ep
#******************************************************************************

#******************************************************************************
# Bollinger Events Finder (events for which the stock price went beyone Bollinger band)
#******************************************************************************
def BollingerEvents(Close_DF, Timestamps, Lookback, Boll_width = -2.0, Benchmark='SPY', Bench_width = 1.5):
    
    print 'Identifying Bollinger Events...'
    Mean_DF = pd.rolling_mean(Close_DF, Lookback) #Dataframe of moving averages of closing prices
    Std_DF = pd.rolling_std(Close_DF, Lookback) #Dataframe of moving standard deviations of closing prices
    Bands_DF = (Close_DF - Mean_DF) / Std_DF #Standardized Dataframe of difference between closing prices and moving averages 
    Events_DF = copy.deepcopy(Close_DF) * np.NAN #initialize the event Dataframe
    Timestamps = Close_DF.index #retrieve time stamps for the event range

    for s_sym in Symbol_List: #loop through symbols in Dataframe
        for i in range(1, len(Timestamps)): #loop through the timeseries
            f_boll_today = Bands_DF[s_sym].ix[Timestamps[i]] #retrieve the current bollinger value (how many standard deviation is the closing price away from moving average)
            f_boll_yest = Bands_DF[s_sym].ix[Timestamps[i - 1]] #retrieve the previous bollinger value (how many standard deviation is the closing price away from moving average)
            f_spy_today = Bands_DF[Benchmark].ix[Timestamps[i]] #look at the current bollinger value for the benchmark
            
            if f_boll_yest >= Boll_width and f_boll_today <= Boll_width and f_spy_today >= Bench_width: #Spotting the event
                Events_DF[s_sym].ix[Timestamps[i]] = 1 #assigning the event to the Event Dataframe
            #end
        #end
    #end
    return Events_DF
#end
#******************************************************************************   

#******************************************************************************
# Orders Routine
#******************************************************************************
def CreateOrders(Events_DF, Sell_Lag = 5, OrderFile='Orders.csv'):
    
    print 'Creating Orders...'
    Out_File = open(OrderFile, 'w') #open the order file for writing
    Out_File.writelines('Year,Month,Day,Ticker,Order,Amount\n') #writes the header

    for symbol in Events_DF.columns: #loop for symbols in the Event Dataframe
        for itime in xrange(len(Timestamps)): #loop through the timeseries
            Buydate = Timestamps[itime] #retrieve the date for buying
            if not np.isnan(Events_DF.get_value(Buydate,symbol)): #check for events presence
                if itime+Sell_Lag >= len(Timestamps): #check whether sell time is out of timeseries
                    Selldate = Timestamps[len(Timestamps) - 1] #assign sell date
                else:
                    Selldate = Timestamps[itime + Sell_Lag] #assign sell date
                #end
                Out_File.writelines(Buydate.strftime('%Y,%m,%d') + ',' + str(symbol) + ',Buy,100\n') #write the Buy order on order file
                Out_File.writelines(Selldate.strftime('%Y,%m,%d') + ',' + str(symbol) + ',Sell,100\n') #write the Sell order on order file
                #end
            #end
        #end
    #end
    Out_File.close()
    return
#end
#******************************************************************************

#******************************************************************************
# Assessing Strategy Potential (Strategy Profiler)
#******************************************************************************    
def AssessStrategy(Events_DF, Data_Dict, Benchmark = 'SPY', Figure_Name='Strategy_Profile.pdf'):
    print 'Assessing Strategy...'
    ep.eventprofiler(Events_DF, Data_Dict, i_lookback=20, i_lookforward=20, s_filename=Figure_Name, b_market_neutral=True, b_errorbars=True, s_market_sym=Benchmark)  
    plt.close()    
    return
#end
#******************************************************************************

#******************************************************************************
# Getting Benchmark Portfolio
#****************************************************************************** 
def BenchmarkPortfolio(StartDate, EndDate, StartCash=100000, Key_List=['close','actual_close'], TimeClose=16, BenchmarkSymbol='$SPX', Provider='Yahoo'):
    Dataobj = da.DataAccess(Provider) #retrieve the data object
    Timestamps = du.getNYSEdays(StartDate, EndDate, dt.timedelta(hours=TimeClose)) #obtain the timestamps
    Symbol_List = [BenchmarkSymbol] #create the list of symbols to retrieve
    List_DF_Data = Dataobj.get_data(Timestamps, Symbol_List, Key_List) #create the list of dataframes
    Data_Dict = dict(zip(Key_List, List_DF_Data)) #build the data dictionary
    Benchmark_DF = Data_Dict['close'].copy() #extract the dataframe for the benchmark
    Relative_Benchmark_DF = Benchmark_DF/Benchmark_DF[BenchmarkSymbol][0] #determine the benchmark relative to the first day of simulation
    BenchmarkPortfolio_DF = Relative_Benchmark_DF * StartCash #evaluate the benchmark portfolio
    return BenchmarkPortfolio_DF
#end
#******************************************************************************

#******************************************************************************
# Creating the Trades Dataframe
#****************************************************************************** 
def TradesDataframe(Data_Dict, Order_DF, Symbol_List, Timestamps, TimeOrder, OrderFee, StartCash=100000):
    Trades_DF = pd.DataFrame(index=list(Timestamps), columns=list(Symbol_List)) #initialize the Dataframe
    Stocks = dict() #initialize the stocks holding portfolio    
    for sym in Symbol_List: #loop thgouht symbols
        Stocks[sym] = 0 #initialize the stock holding at time 0
        Trades_DF[sym][Timestamps[0]] = 0 #initialize the trades at time 0
    #end
    Cash = StartCash #initialize the current cash
    Trades_DF['CASH'][Timestamps[0]] = StartCash #initialize the portfolio in the Trades Dataframe    
    for index, row in Order_DF.iterrows(): #loop through orders list
        Date = dt.datetime(row['Year'], row['Month'], row['Day'], TimeOrder ) #retrieve the date of the order
        Sym = row['Ticker'] #retrieve the stock symbol of the order
        Stock_Value = Data_Dict['close'][Sym][Date] #retrieve the stock price at the time of the order
        Stock_Amount = row['Amount'] #retrieve the stock amount
        if row['Order'] == 'Buy': #if Buy order
            Cash = Cash - (Stock_Value * Stock_Amount) - OrderFee #calculate the new cash balance
            Trades_DF['CASH'][Date] = Cash #update the balance in the Trades Dataframe
            Stocks[Sym] = Stocks[Sym] + Stock_Amount #calculate the number of stock held
            Trades_DF[Sym][Date] = Stocks[Sym] #update the number of stocks held in theTrades Dataframe
        else: #if Sell order
            Cash = Cash + (Stock_Value * Stock_Amount) - OrderFee #calculate the new cash balance
            Trades_DF['CASH'][Date] = Cash #update the balance in the Trades Dataframe
            Stocks[Sym] = Stocks[Sym] - Stock_Amount #calculate the number of stock held
            Trades_DF[Sym][Date] = Stocks[Sym] #update the number of stocks held in theTrades Dataframe
        #end            
    #end
    Trades_DF = Trades_DF.fillna(method = 'pad') #handling missing data (padding the remaining rows (in between orders nothing changes in the portfolio))
    return Trades_DF
#end
#******************************************************************************
    
#******************************************************************************
# Creating the Portfolio Value Dataframe
#****************************************************************************** 
def PortfolioValueDF(Data_Dict, Trades_DF, Symbol_List, Timestamps):
    Value_DF = pd.DataFrame(index=list(Timestamps), columns=['Value']) #create the Dataframe with the Portfolio value
    Value_DF = Value_DF.fillna(0) #fill the Dataframe (initialize with 0s)

    for day in Timestamps: #loop through timestamps
        Value = 0 #initialize the portfolio value
        for sym in Symbol_List: #loop thoguh symbols (stocks)
            if sym == 'CASH': #for Cash balance
                Value = Value + Trades_DF[sym][day] #add the cash balance to the portfolio value
            else: #for all other symbols (stocks)
                Value = Value + Trades_DF[sym][day] * Data_Dict['close'][sym][day] #update the portfolio value by looking at the value of the stocks held
            #end
        #end
        Value_DF['Value'][day] = Value #update the daily value of the portfolio
    #end
    return Value_DF
#end
#******************************************************************************
    
#******************************************************************************
# Portfolio Simulator
#******************************************************************************
def PortfolioSimulator(Data_Dict, StartDate, EndDate, OrderFile='Orders.csv', TimeOrder=16, ValueFile='Portfolio.csv', StartCash=100000, OrderFee=0):

    print 'Simulating Portfolio...'
    
    #------------------------------------------------------
    # Reading from Order file
    #------------------------------------------------------        
    Order_DF = pd_par.read_csv(OrderFile, header=0).sort(['Year','Month','Day']) #read the Order Dataframe and sorting the trades DF by increasing date    
    Symbol_List = list(set(Order_DF['Ticker'].values)) #getting the Symbols from the .csv file
    Symbol_List.append('CASH') #adding CASH to list of Symbols and creating the trades table
    
    #------------------------------------------------------
    # Creating a Trades Dataframe containing the list of all trades and portfolio holdings
    #------------------------------------------------------    
    Trades_DF = TradesDataframe(Data_Dict, Order_DF, Symbol_List, Timestamps, TimeOrder, OrderFee, StartCash=100000) #call routine to create Trades Dataframe
    #------------------------------------------------------
    
    #------------------------------------------------------
    # Creating a Dataframe with the portfolio value
    #------------------------------------------------------     
    Value_DF = PortfolioValueDF(Data_Dict, Trades_DF, Symbol_List, Timestamps)  #call routine to create Portfolio Value Dataframe
    #------------------------------------------------------
    
    #------------------------------------------------------
    # Retrieve the Benchmark Portfolio
    #------------------------------------------------------
    BenchmarkSymbol = '$SPX' #set the benchmark symbol ($SPX is the ETF for the S&P500)
    BenchmarkPortfolio_DF = BenchmarkPortfolio(StartDate,EndDate,StartCash,Key_List =['close','actual_close'], TimeClose=16, BenchmarkSymbol=BenchmarkSymbol) #retrieve the benchmark portfolio
    #------------------------------------------------------
    
    #------------------------------------------------------
    # Save the output
    #------------------------------------------------------
    Out_File = open(ValueFile, 'w') #open the output file for writing
    Out_File.writelines('Year,Month,Day,Portfolio Value,Benchmark Value\n') #writes the header
    for index, row in Value_DF.iterrows(): #loop through portfolio Values Dataframe
        Out_File.writelines(str(index.strftime('%Y,%m,%d')) + ',' + str(row['Value'].round()) + ',' + str(BenchmarkPortfolio_DF[BenchmarkSymbol][index].round()) + '\n' ) #write lines on output file
    #end
    Out_File.close() #close the output file
    #------------------------------------------------------
    
    return Trades_DF, Value_DF #return the Trades and portfolio Value Dataframes
    
#end
#******************************************************************************

#******************************************************************************
# Plot the Portfolio and Corresponding Benchmark
#******************************************************************************
def PlotPortfolio(Value_DF, Benchmark_DF, FigureName='Simulated Portfolio.pdf', Fig_counter=1):
    plt.figure(Fig_counter)
    Value_DF['Value'].plot(title='Portfolio Simulation', linewidth=2, label='Event Strategy', legend=True)
    Benchmark_DF['$SPX'].plot(title='Portfolio Simulation', linewidth=2, label='S&P500 Benchmark', legend=True)
    plt.ylabel('Portfolio Value ($)')
    plt.savefig(FigureName)
    plt.close(Fig_counter)
    return
#end
#******************************************************************************

#******************************************************************************
# Main Program
#******************************************************************************

if __name__ == '__main__': #execute the main program if the module is called directly and not as part of an import

    #-----------------------
    # Initialization
    #-----------------------
    StartDate = dt.datetime(2008, 1, 1) #define the start date
    EndDate = dt.datetime(2009, 12, 31) #define the end date
    Lookback = 20 #lookback 20 days in strategy evaluation
    TimeClose = 16 #define the time of the day for stock value retrieval
    Provider = 'Yahoo'
    SymbolsCategory = 'sp5002012'
    BenchmarkSymb = 'SPY'
    Key_List = ['close','actual_close']
    #-----------------------
    
    #-----------------------
    # Retrieve Data
    #-----------------------
    print 'Retrieving Data...'
    Dataobj = da.DataAccess(Provider) #retrieve the data object from QSTK module
    Symbol_List = Dataobj.get_symbols_from_list(SymbolsCategory) #retrieve the list of symbols
    Symbol_List.append(BenchmarkSymb) #append the benchmark symbol to the list of symbols
    Timestamps = du.getNYSEdays(StartDate, EndDate, dt.timedelta(hours=TimeClose)) #obtain the timestams from QSTK module
    List_DF_data = Dataobj.get_data(Timestamps, Symbol_List, Key_List) #build a list of Dataframes for each Key and each Symbol 
    Data_Dict = dict(zip(Key_List, List_DF_data)) #turn the Dataframe into a dictionary, using the tuple of Key and Dataframe for that Key
    #-----------------------
    
    #-----------------------
    # Handling Missing Data
    #-----------------------  
    for s_key in Key_List: #loop through keys
        Data_Dict[s_key] = Data_Dict[s_key].fillna(method = 'ffill') #forward fill NaN
        Data_Dict[s_key] = Data_Dict[s_key].fillna(method = 'bfill') #backward fill NaN
        Data_Dict[s_key] = Data_Dict[s_key].fillna(1.0) #fill remaining NaN
    #end
    Close_DF = Data_Dict['close'] #retrieve the Dataframe for closing prices of all symbols
    #-----------------------

    #-----------------------
    # Creating Output
    #-----------------------  
    Events_DF = BollingerEvents(Close_DF,Timestamps,Lookback) #Find the Bollinger Events occurring in the study period
    AssessStrategy(Events_DF, Data_Dict, Benchmark = 'SPY', Figure_Name='Strategy_Profile.pdf') #assessing strategy (plotting the profiler)
    CreateOrders(Events_DF, Sell_Lag = 5, OrderFile='Orders.csv') #create the order file for the strategy
    Trades_DF, Values_DF = PortfolioSimulator(Data_Dict, StartDate, EndDate, OrderFile='Orders.csv', ValueFile='Portfolio.csv', StartCash=100000, OrderFee=0) #Simulating the Portfolio   
    Benchmark_DF = BenchmarkPortfolio(StartDate,EndDate,StartCash=100000,Key_List=['close','actual_close'], TimeClose=16, BenchmarkSymbol='$SPX') #retrieve the benchmark portfolio   
    PlotPortfolio(Values_DF,Benchmark_DF, FigureName='Simulated Portfolio.pdf', Fig_counter=2) #plot portfolio and benchmark and save the figure
    #-----------------------
    
#end
#******************************************************************************

