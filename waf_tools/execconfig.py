import waflib.TaskGen
import waflib.Task as Task
from waflib import Utils

@waflib.TaskGen.feature("*")
@waflib.TaskGen.before_method('process_source')
def process_execrule(self):
  
  if not getattr(self,'execrule',None):
    return
  self.meths.remove('process_source')
  name=str(getattr(self,'name',None)or self.target or self.execrule)
  cls=Task.task_factory(name,self.execrule,getattr(self,'vars',[]),shell=getattr(self,'shell',True),color=getattr(self,'color','BLUE'))
  tsk=self.create_task(name)
  if getattr(self,'target',None):
    if isinstance(self.target,str):
      self.target=self.target.split()
    if not isinstance(self.target,list):
      self.target=[self.target]
    for x in self.target:
      if isinstance(x,str):
        tsk.outputs.append(self.path.find_or_declare(x))
      else:
        x.parent.mkdir()
        tsk.outputs.append(x)
    if getattr(self,'install_path',None):
      self.bld.install_files(self.install_path,tsk.outputs,chmod=Utils.O755)
  if getattr(self,'source',None):
    tsk.inputs=self.to_nodes(self.source)
    self.source=[]
  if getattr(self,'scan',None):
    cls.scan=self.scan
  if getattr(self,'cwd',None):
    tsk.cwd=self.cwd
  if getattr(self,'update_outputs',None)or getattr(self,'on_results',None):
    Task.update_outputs(cls)
  if getattr(self,'always',None):
    Task.always_run(cls)
  for x in['after','before','ext_in','ext_out']:
    setattr(cls,x,getattr(self,x,[]))

def mkstr(pat,li):
  return " ".join([pat%l for l in li])
def bld_link_line(bld,name):
  if name:
    return mkstr(bld.env.LIBPATH_ST,getattr(bld.env,"LIBPATH_"+name)) + " " +   mkstr(bld.env.RPATH_ST,getattr(bld.env,"LIBPATH_"+name)) + " " + mkstr(bld.env.LIB_ST,getattr(bld.env,"LIB_"+name))
  return ""
def bld_inc_line(bld,name):
  if name:
    return mkstr(bld.env.CPPPATH_ST,getattr(bld.env,"INCLUDES_"+name))
  return ""
def bld_def_line(bld,name):
  if name:
    return mkstr(bld.env.DEFINES_ST,getattr(bld.env,"DEFINES_"+name))
  return ""

conftemplate = """#! /usr/bin/env python
# don't do much for now
print "%s"
"""
def createconfig(*confline):
  def doconf(task):
    tgt = task.outputs[0].abspath() 
    f=open(tgt,"w")
    print >>f,conftemplate%" ".join(confline),
    f.close()  
  return doconf

