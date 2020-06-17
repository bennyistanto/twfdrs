function FFMCcalc(ffmc0,prec,hum,wind,temp){
  var mo = ffmc0.expression('(147.2*(101.0 - b(0)))/(59.5 + b(0))') //*Eq. 1* 
  mo = img.expression('prec>0.5? mo>150? eq3b : eq3a : mo',{
    'prec':prec.select(0),
    'mo':mo.select(0),
    'eq3b':mo.expression('(mo + 42.5 * rf * exp(-100.0/(251.0-mo))*(1.0 - exp(-6.93/rf))) + (0.0015*pow((mo - 150.0),2))*sqrt(rf)',{
      'mo':mo.select(0),
      'exp':Math.exp(),
      'pow':Math.pow(),
      'sqrt':Math.sqrt(),
      'rf':prec.subtract(0.5)}),
    'eq3a':mo.expression('mo + 42.5*rf*exp(-100.0/(251.0-mo))*(1.0 - exp(-6.93/rf))',{
      'mo':mo.select(0),
      'exp':Math.exp(),
      'sqrt':Math.sqrt(),
      'rf':prec.subtract(0.5)})})
  var ed = img.expression('0.942*(pow(hum,0.679)) + (11.0*exp((hum-100.0)/10.0))+ 0.18*(21.1-temp)*(1.0 - 1.0/exp(0.1150 * hum))',{
    'pow':Math.pow(),
    'exp':Math.exp(),
    'hum':hum.select(0),
    'temp':temp.select(0)}) //*Eq. 4*#
  var m = img.expression('mo<ed? mo<=ew? ew - (ew - mo)/(pow(10.0,kw1)):mo : mo==ed? mo: ed + (mo-ed)/(pow(10.0,kw2))',{
    'pow':Math.pow(),
    'mo':mo.select(0),
    'ed':ed.select(0),
    'ew':img.expression('0.618*(pow(hum,0.753)) + (10.0*exp((hum-100.0)/10.0))+ 0.18*(21.1-temp)*(1.0 - 1.0/exp(0.115 * hum))',{
      'pow':Math.pow(),
      'exp':Math.exp(),
      'hum':hum.select(0),
      'temp':temp.select(0)}),
    'kw1':temp.expression('kl1 * (0.581 * exp(0.0365 * temp))',{
      'exp':Math.exp(),
      'temp':temp.select(0),
      'kl1':wind.expression('0.424*(1.0-pow(((100.0-hum)/100.0),1.7))+(0.0694*sqrt(wind))*(1.0 - pow(((100.0 - hum)/100.0),8))',{
        'pow':Math.pow(),
        'sqrt':Math.sqrt(),
        'hum':hum.select(0),
        'wind':wind.select(0)})}),
    'kw2':img.expression('kl2 * (0.581 * exp(0.0365 * temp))',{
      'exp':Math.exp(),
      'temp':temp.select(0),
      'kl2':wind.expression('0.424*(1.0-pow((hum/100.0),1.7))+(0.0694*sqrt(wind))*(1.0-pow((hum/100.0),8))',{
        'pow':Math.pow(),
        'sqrt':Math.sqrt(),
        'hum':hum.select(0),
        'wind':wind.select(0)})})})
  var ffmc = m.expression('(59.5 * (250.0 - b(0))) / (147.2 + b(0))')
  ffmc = ffmc.expression('b(0)>101.0? 101.0 : b(0)<=0.0? 0.0 : b(0)')
  return ffmc }

function DMCcalc(dmc0,mth,temp,prec,hum){
  var el = [6.5,7.5,9.0,12.8,13.9,13.9,12.4,10.9,9.4,8.0,7.0,6.0]
  var t = temp.expression('b(0)<-1.1? -1.1:b(0)')
  var rk = img.expression('1.894*(temp+1.1) * (100.0-hum) * (el*0.0001)',{
    'temp':t.select(0),
    'hum':hum.select(0),
    'el':el[mth-1]})
  var b = dmc0.expression('b(0)<=33.0? eq13a : b(0)<=65? eq13b : eq13c',{
    'eq13a':dmc0.expression('100.0 /(0.5 + 0.3*b(0))'),
    'eq13b':dmc0.expression('14.0 - 1.3*log(b(0))',{'log':Math.log()}),
    'eq13c':dmc0.expression('6.2 * log(b(0)) - 17.2',{'log':Math.log()})})
  var pr = img.expression('prec>1.5? 43.43*(5.6348 - log(wmr-20.0)) : dmc0',{
    'log':Math.log(0),
    'prec':prec.select(0),
    'dmc0':dmc0.select(0),
    'wmr':img.expression('wmi + (1000*rw) / (48.77+b*rw)',{
      'wmi':img.expression('20.0 + 280.0/exp(0.023*dmc0)',{
        'exp':Math.exp(),
        'dmc0':dmc0.select(0)}),
      'rw':prec.expression('0.92*b(0)-1.27'),
      'b':b.select(0)})})
  pr = pr.expression('b(0)<0?0:b(0)')
  var dmc = img.expression('pr+rk',{
    'pr':pr.select(0),
    'rk':rk.select(0)})
  dmc = dmc.expression('b(0)<=1.0?1.0:b(0)')
  return dmc}

function DCcalc(dc0,mth,temp,prec){
  var fl = [-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6]
  var t = temp.expression('b(0)<-2.8?-2.8:b(0)')
  var pe = t.expression('(0.36*(t+2.8)+fl)/2',{
    't':t.select(0),
    'fl':fl[mth-1]})
  pe = pe.expression('b(0)<=0.0?0.0:b(0)')
  var dc = img.expression('prec>2.8? dr>0.0? dr+pe: pe : dc0+pe',{
    'prec':prec.select(0),
    'pe':pe.select(0),
    'dc0':dc0.select(0),
    'dr':dc0.expression('dc0 - 400.0*log( 1.0+((3.937*rw)/smi))',{
      'log':Math.log(),
      'dc0':dc0.select(0),
      'smi':dc0.expression('800.0*exp(-b(0)/400.0)',{'exp':Math.exp()}),
      'rw':prec.expression('0.83*b(0) - 1.27')
    })})
  return dc }
  
function ISIcalc(ffmc,wind){
  var mo = ffmc.expression('147.2*(101.0-b(0))/(59.5+b(0))') //#*Eq. 1*#
  var ff = mo.expression('19.115*exp(b(0)*-0.1386)*(1.0+(pow(b(0),5.31))/49300000.0)',{
    'exp':Math.exp(),
    'pow':Math.pow()}) //*Eq. 25*#
  var isi = ff.expression('b(0)*exp(0.05039*wind)',{
    'exp':Math.exp(),
    'wind':wind.select(0)}) //*Eq. 26*#
  return isi}

function BUIcalc(dmc,dc){
  var bui  = img.expression('dmc<=0.4*dc? (0.8*dc*dmc)/(dmc+0.4*dc)'+
  ': dmc-(1.0-0.8*dc/(dmc+0.4*dc))*(0.92+pow((0.0114*dmc),1.7))',{
    'pow':Math.pow(),
    'dmc':dmc.select(0),
    'dc':dc.select(0)})
  bui = bui.expression('b(0)<0.0?0.0:b(0)')
  return bui}
  
function FWIcalc(isi,bui){
  var bb = bui.expression('b(0)<=80.0? 0.1*isi*(0.626*pow(bui,0.809)+2.0)'+
  ': 0.1*isi*(1000.0/(25.0 + 108.64/exp(0.023*bui)))',{
    'pow':Math.pow(),
    'exp':Math.exp(),
    'bui':bui.select(0),
    'isi':isi.select(0)})
  var fwi = bb.expression('b(0)<=1.0?b(0):exp(2.72 * pow((0.434*log(b(0))),0.647))',{
    'exp':Math.exp(),
    'pow':Math.pow()})
  return fwi}
  
function main(ffmc0,dmc0,dc0,prec,hum,wind,temp,mth){
  var ffmc = FFMCcalc(ffmc0,prec,hum,wind,temp)
  var dmc = DMCcalc(dmc0,mth,temp,prec,hum)
  var dc = DCcalc(dc0,mth,temp,prec)
  var isi = ISIcalc(ffmc,wind)
  var bui = BUIcalc(dmc,dc)
  var fwi = FWIcalc(isi,bui)
  return print(ffmc,dmc,dc,isi,bui,fwi)
}

var ffmc0 = ee.Image.constant(85)
var dmc0 = ee.Image.constant(6)
var dc0 = ee.Image.constant(15)
var prec = ee.Image.random()
var hum = ee.Image.constant(42)
var wind = ee.Image.constant(25)
var temp = ee.Image.constant(17)
var mth = 4
var img = ee.Image.constant(0)

main(ffmc0,dmc0,dc0,prec,hum,wind,temp,mth)
