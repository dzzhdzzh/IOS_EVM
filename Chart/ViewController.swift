//
//  ViewController.swift
//  Chart
//
//  Created by Bruno Sielly Jales Costa on 8/2/16.
//  Copyright © 2016 Hackathon. All rights reserved.
//

import UIKit
import Charts
import CoreBluetooth
import AVFoundation

let WINDOW_SIZE = 20;

class ViewController: UIViewController, ChartViewDelegate, CBCentralManagerDelegate, CBPeripheralDelegate, CvVideoCameraDelegate{

    @IBOutlet var lineChartView: LineChartView!;
    @IBOutlet var uiHeart: UIImageView!;
    @IBOutlet var hbLabel: UILabel!;
    //@IBOutlet var previewView: UIView!
    @IBOutlet var uiHeart2: UIImageView!;
    @IBOutlet var hb2Label: UILabel!;
    @IBOutlet var btOn: UISwitch!;
    @IBOutlet var previewView: UIImageView!;
    
    var centralManager: CBCentralManager?;
    var peripherals: Array<CBPeripheral> = Array<CBPeripheral>();
    var BTHeartMonitor: CBPeripheral?;
    var BTHeartRate: Double = 0;
    var EstimatedHeartRate: Double = 0;
    var captureSession: AVCaptureSession?;
    var stillImageOutput: AVCaptureStillImageOutput?;
    var previewLayer: AVCaptureVideoPreviewLayer?;
    var currentView: UIView?;
    var myCamera : CvVideoCamera!;
    var videoCameraWrapper : CvVideoCameraWrapper!;
    var meanHR: Double = 0;
    var countHR: Int8 = 0;
    
    
    override func viewDidLoad() {
        super.viewDidLoad();
        centralManager = CBCentralManager(delegate: self, queue: nil);
        self.lineChartView.delegate = self;
        self.lineChartView.descriptionTextColor = UIColor.whiteColor();
        self.lineChartView.gridBackgroundColor = UIColor.darkGrayColor();
        self.lineChartView.noDataText = "No data provided";
        setChartData();
        self.currentView = self.previewView;
    }

    override func didReceiveMemoryWarning() {
        super.didReceiveMemoryWarning();
        // Dispose of any resources that can be recreated.
    }
    
    @IBAction func btnOriginal(sender: AnyObject) {
        //if (self.currentView != self.previewView){
            self.view.bringSubviewToFront(self.previewView);
            self.currentView = self.previewView;
            self.videoCameraWrapper.setOriginal(1);
        //}
        
    }
    
    @IBAction func btnHistory(sender: AnyObject) {
        //if (self.currentView != self.lineChartView){
            self.view.bringSubviewToFront(self.lineChartView);
            self.currentView = self.lineChartView;
        //}
        print("SWIFT ", self.videoCameraWrapper.getHR());
    }
    
    @IBAction func btnMagnified(sender: AnyObject) {
        self.view.bringSubviewToFront(self.previewView);
        self.currentView = self.previewView;
        self.videoCameraWrapper.setOriginal(0);
    }
    
    @IBAction func btSelected(sender: UISwitch) {
        if (!sender.on){
            self.uiHeart2.alpha = 0.0;
            self.hb2Label.alpha = 0.0;
            self.lineChartView.data?.dataSets[1].visible = false;
        }
        else{
            self.lineChartView.data?.dataSets[1].visible = true;
        }
    }
        
    override func viewWillAppear(animated: Bool) {
        super.viewWillAppear(animated);
        self.videoCameraWrapper = CvVideoCameraWrapper(controller:self, andImageView: previewView);
        self.videoCameraWrapper.setOriginal(1);
    }
    
    func updateChart() {
        let data: ChartData = self.lineChartView.data!;
        let dataSet1: IChartDataSet = data.getDataSetByIndex(0);
        let dataSet2: IChartDataSet = data.getDataSetByIndex(1);
        dataSet1.addEntry(ChartDataEntry(value: self.EstimatedHeartRate, xIndex: dataSet1.entryCount));
        dataSet2.addEntry(ChartDataEntry(value: self.BTHeartRate, xIndex: dataSet2.entryCount));
        data.addXValue(String(data.xValCount + 1));
        data.dataSets[0] = dataSet1;
        data.dataSets[1] = dataSet2;
        self.lineChartView.data = data;
        self.lineChartView.setVisibleXRangeMaximum(CGFloat(WINDOW_SIZE - 1));
        self.lineChartView.moveViewToX(CGFloat(dataSet1.entryCount));
    }
    
    func animateHeart(){
        //self.EstimatedHeartRate = Double(arc4random_uniform(10) + 70);
        self.EstimatedHeartRate = self.videoCameraWrapper.getHR();
        self.uiHeart.alpha = 1.0;
        self.hbLabel.alpha = 1.0;
        self.hbLabel.text = String(Int16(self.EstimatedHeartRate));
        if (self.EstimatedHeartRate > 0){
            UIView.animateWithDuration(60.0/self.EstimatedHeartRate, delay: 0, options: UIViewAnimationOptions.CurveEaseOut, animations:{ self.uiHeart.alpha = 0.3; }, completion: nil);
        }
        else{
            //UIView.animateWithDuration(5, delay: 0, options: UIViewAnimationOptions.CurveEaseOut, animations:{ self.uiHeart.alpha = 0.3; }, completion: nil);
        }
    }
    
    func setChartData() {
        
        self.lineChartView.descriptionText = "";
        self.lineChartView.leftYAxisRenderer.yAxis?.labelTextColor = UIColor.whiteColor();
        self.lineChartView.legend.textColor = UIColor.whiteColor();
        self.lineChartView.setViewPortOffsets(left: 30, top: 40, right: 10, bottom: 10)
        self.lineChartView.legend.enabled = false;
        
        let yVals1 : [ChartDataEntry] = [ChartDataEntry]();
        let set1: LineChartDataSet = LineChartDataSet(yVals: yVals1, label: "Estimated Heart Rate");
        
        set1.axisDependency = .Left;
        set1.setColor(UIColor.redColor().colorWithAlphaComponent(1));
        set1.lineWidth = 2.0;
        set1.circleRadius = 0.0; // the radius of the node circle
        set1.fillAlpha = 65 / 255.0;
        set1.fillColor = UIColor.redColor();
        set1.highlightEnabled = false;
        set1.drawCircleHoleEnabled = false;
        set1.drawCirclesEnabled = false;
        set1.drawValuesEnabled = false;
        
        let yVals2 : [ChartDataEntry] = [ChartDataEntry]();
        let set2: LineChartDataSet = LineChartDataSet(yVals: yVals2, label: "BT Heart Monitor");
        
        set2.axisDependency = .Left;
        set2.setColor(UIColor(red: 0.082, green: 0.396, blue: 0.953, alpha: 1));
        set2.lineWidth = 2.0;
        set2.circleRadius = 0.0; // the radius of the node circle
        set2.fillAlpha = 65 / 255.0;
        set2.fillColor = UIColor.redColor();
        set2.highlightEnabled = false;
        set2.drawCircleHoleEnabled = false;
        set2.drawCirclesEnabled = false;
        set2.drawValuesEnabled = false;
        
        //3 - create an array to store our LineChartDataSets
        var dataSets : [LineChartDataSet] = [LineChartDataSet]();
        let months : [Int] = [Int]();
        dataSets.append(set1);
        dataSets.append(set2);
        set2.visible = false;
        
        //4 - pass our months in for our x-axis label value along with our dataSets
        let data: LineChartData = LineChartData(xVals: months, dataSets: dataSets)
        data.setValueTextColor(UIColor.whiteColor())
        
        //5 - finally set our data
        self.lineChartView.data = data
        self.lineChartView.setVisibleXRangeMinimum(0);
        self.lineChartView.setVisibleXRangeMaximum(20);
        self.lineChartView.leftAxis.axisMaxValue = 180;
        self.lineChartView.leftAxis.axisMinValue = 30;
        self.lineChartView.rightAxis.axisMaxValue = 180;
        self.lineChartView.rightAxis.axisMinValue = 30;
        
        _ = NSTimer.scheduledTimerWithTimeInterval(1.0, target: self, selector: #selector(ViewController.updateChart), userInfo: nil, repeats: true);
        _ = NSTimer.scheduledTimerWithTimeInterval(1.0, target: self, selector: #selector(ViewController.animateHeart), userInfo: nil, repeats: true);
    }
    
    func centralManagerDidUpdateState(central: CBCentralManager) {
        print("centralManagerDidUpdateState")
        if (central.state == CBCentralManagerState.PoweredOff) {
            print("CoreBluetooth BLE hardware is powered off")
        }
        else if (central.state == CBCentralManagerState.PoweredOn) {
            print("CoreBluetooth BLE hardware is powered on and ready")
            centralManager!.scanForPeripheralsWithServices(nil, options: nil)
        }
        else if (central.state == CBCentralManagerState.Unauthorized) {
            print("CoreBluetooth BLE state is unauthorized")
        }
        else if (central.state == CBCentralManagerState.Unknown) {
            print("CoreBluetooth BLE state is unknown")
        }
        else if (central.state == CBCentralManagerState.Unsupported) {
            print("CoreBluetooth BLE hardware is unsupported on this platform")
        }
    }
    
    
    func centralManager(central: CBCentralManager, didDiscoverPeripheral peripheral: CBPeripheral, advertisementData: [String : AnyObject], RSSI: NSNumber){
        let deviceName = "Polar H7 BA784219";
        let nameOfDeviceFound = (advertisementData as NSDictionary).objectForKey(CBAdvertisementDataLocalNameKey) as? NSString;
        if (nameOfDeviceFound == deviceName){
            print("found " + String(nameOfDeviceFound));
            self.BTHeartMonitor = peripheral;
            self.BTHeartMonitor!.delegate = self;
            centralManager?.connectPeripheral(peripheral, options: nil);
            centralManager?.stopScan();
            print(peripheral);
        }
        else{
            print("-> not found");
        }
    }
    
    func centralManager(central: CBCentralManager, didFailToConnectPeripheral peripheral: CBPeripheral, error: NSError?) {
        print("-> failed to connect");
    }
    
    func centralManager(central: CBCentralManager, didConnectPeripheral peripheral: CBPeripheral) {
        print("-> connected to device");
        peripheral.discoverServices(nil);
        
        //UIView.animateWithDuration(1.5, delay: 3.0, options: UIViewAnimationOptions.CurveEaseOut, animations: {
            //self.uiHeart.alpha = 1.0;
            //}, completion: nil)
    }
    
    func centralManager(central: CBCentralManager, didDisconnectPeripheral peripheral: CBPeripheral, error: NSError?) {
        print("-> disconnected");
    }
    
    func peripheral(peripheral: CBPeripheral, didDiscoverServices error: NSError?) {
        for service in peripheral.services! {
            let thisService = service as CBService;
            print(thisService.UUID);
            if String(thisService.UUID) == String("Heart Rate") {
                print("-> service discovered ");
                self.BTHeartMonitor!.discoverCharacteristics(nil, forService: thisService);
            }
        }
    }
    
    func peripheral(peripheral: CBPeripheral, didDiscoverCharacteristicsForService service: CBService, error: NSError?) {
        for charateristic in service.characteristics! {
            let thisCharacteristic = charateristic as CBCharacteristic;
            if String(thisCharacteristic.UUID) == String("2A37") {
                peripheral.setNotifyValue(true, forCharacteristic: thisCharacteristic)
                peripheral.readValueForCharacteristic(thisCharacteristic);
                self.BTHeartMonitor = peripheral;
                print("-> setting notification");
            }
        
        }
    }
    
    func peripheral(peripheral: CBPeripheral, didUpdateValueForCharacteristic characteristic: CBCharacteristic, error: NSError?) {
        if String(characteristic.UUID) == String("2A37") {
            let dataBytes = characteristic.value;
            
            //let dataLength = dataBytes!.length;
            //var dataArray = [Int16](count: dataLength, repeatedValue: 0);
            //dataBytes!.getBytes(&dataArray, length: dataLength * sizeof(Int16));
            
            var buffer = [UInt8](count: dataBytes!.length, repeatedValue: 0x00);
            dataBytes!.getBytes(&buffer, length: buffer.count)
            var bpm:UInt16?
            if (buffer.count >= 2){
                if (buffer[0] & 0x01 == 0){
                    bpm = UInt16(buffer[1]);
                }else {
                    bpm = UInt16(buffer[1]) << 8
                    bpm =  bpm! | UInt16(buffer[2])
                }
            }
            self.BTHeartRate = Double(bpm!);
            
            if (btOn.on){
                self.uiHeart2.alpha = 1.0;
                self.hb2Label.alpha = 1.0;
                self.hb2Label.text = String(Int16(self.BTHeartRate));
                if (self.BTHeartRate > 0){
                    UIView.animateWithDuration(60.0/self.BTHeartRate, delay: 0, options: UIViewAnimationOptions.CurveEaseOut, animations:{ self.uiHeart2.alpha = 0.3; }, completion: nil);
                }
                else{
                    UIView.animateWithDuration(1, delay: 0, options: UIViewAnimationOptions.CurveEaseOut, animations:{ self.uiHeart2.alpha = 0.3; }, completion: nil);
                }
            }
            
        }
    }

}