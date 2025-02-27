# Copyright (c) 2013, 2014 All Right Reserved, Hatfull Lab, University of Pittsburgh
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
# PARTICULAR PURPOSE.  USE AT YOUR OWN RISK.
#
# Marissa Pacey
# April 4, 2014
# Starterate Window

import MySQLdb
from gi.repository import Gtk, Gdk
import threading
import time
from .uiDialogs import StarteratorFinishedDialog
from . import starterate
from .utils import StarteratorError
from . import database
from . import phamgene

class StarteratorEnterInformation(Gtk.Dialog):
    def show_starterate_button(self):
        self.starterate_button = Gtk.Button("Starterate")
        self.starterate_button.connect('clicked', self.starterate)
        self.box.pack_start(self.starterate_button, False, False, 0)

    def show_phage_entry(self):
        phage_label = Gtk.Label("Phage:")
        phage_list_entry = Gtk.Entry()
        phage_list_entry.connect("changed", self.on_entry_changed, 'phage')
        hbox = Gtk.Box(spacing= 6)
        hbox.pack_start(phage_label, False, False, 0)
        hbox.pack_start(phage_list_entry, False, False, 0)
        self.box.pack_start(hbox, False, False, 0)

    def show_gene_no_entry(self):
        hbox = Gtk.Box(spacing= 6)
        gene_no_label = Gtk.Label("Gene Number:")
        self.gene_no_entry = Gtk.Entry()
        self.gene_no_entry.connect("changed", self.on_entry_changed, 'gene_no')
        hbox.pack_start(gene_no_label, False, False, 0)
        hbox.pack_start(self.gene_no_entry, False, False, 0)
        self.box.pack_start(hbox, False, False, 0 )
       

    def show_fasta_entry(self):
        hbox = Gtk.Box(spacing= 6)
        label = Gtk.Label("Phage Fasta File")
        unphamed_fasta_button = Gtk.Button('Select')
        fasta_entry = Gtk.Entry()
        fasta_entry.connect('changed', self.on_entry_changed, 'fasta')
        unphamed_fasta_button.connect('clicked', self.on_file_clicked, ['fasta', fasta_entry])
        hbox.pack_start(label, False, False, 0)
        hbox.pack_start(fasta_entry, False, False, 0 )
        hbox.pack_start(unphamed_fasta_button, False, False, 0 )
        self.box.pack_start(hbox, False, False, 0)

    def show_profile_entry(self):
        hbox = Gtk.Box(spacing= 6)
        label = Gtk.Label("Phage Profile File (optional)")
        unphamed_profile_button = Gtk.Button('Select')
        profile_entry = Gtk.Entry()
        profile_entry.connect('changed', self.on_entry_changed, 'profile')
        unphamed_profile_button.connect('clicked', self.on_file_clicked, ['profile', profile_entry])
        hbox.pack_start(label, False, False, 0)
        hbox.pack_start(profile_entry, False, False, 0 )
        hbox.pack_start(unphamed_profile_button, False, False, 0 )
        self.box.pack_start(hbox, False, False, 0 )


    def show_gene_info_entry(self):
        hbox = Gtk.Box(spacing= 6)
        start_label = Gtk.Label("Start")
        self.given_start_entry = Gtk.Entry()
        self.given_start_entry.connect("changed", self.on_entry_changed, 'start')
        hbox.pack_start(start_label, False, False, 0)
        hbox.pack_start(self.given_start_entry, False, False, 0)
        self.box.pack_start(hbox, False, False, 0)
        
        hbox = Gtk.Box(spacing= 6)
        stop_label = Gtk.Label("Stop")
        self.given_stop_entry = Gtk.Entry()
        self.given_stop_entry.connect("changed", self.on_entry_changed, 'stop')
        hbox.pack_start(stop_label, False, False, 0)
        hbox.pack_start(self.given_stop_entry, False, False, 0)
        self.box.pack_start(hbox, False, False, 0)
                
        hbox = Gtk.Box(spacing= 6)
        orientation_label = Gtk.Label('Orientation')
        self.orientation_F = Gtk.CheckButton(label="Forward")
        self.orientation_R = Gtk.CheckButton(label="Reverse")
        self.orientation_F.connect("toggled", self.on_orientation_toggled, 'F')
        self.orientation_R.connect("toggled", self.on_orientation_toggled, 'R')
        hbox.pack_start(orientation_label, False, False, 0)
        hbox.pack_start(self.orientation_F, False, False, 0)
        hbox.pack_start(self.orientation_R, False, False, 0)
        self.box.pack_start(hbox, False, False, 0)

    def show_one_pham_entry(self):
        pham_label = Gtk.Label("Pham Number:")
        pham_entry = Gtk.Entry()
        pham_entry.connect("changed", self.on_entry_changed, 'pham')
        hbox = Gtk.Box(spacing= 6)
        hbox.pack_start(pham_label, False, False, 0)
        hbox.pack_start(pham_entry, False, False, 0)
        self.box.pack_start(hbox, False, False, 0)

    def phamerated_all(self):
        self.info['phamerated'] = True
        self.info['all'] = True
        self.show_phage_entry()
        self.show_starterate_button()
        self.show_all()

    def unphamerated_all(self):
        self.info['phamerated'] = False
        self.info['all'] = True
        self.show_phage_entry()
        self.show_fasta_entry()
        self.show_profile_entry()
        self.show_starterate_button()
        self.show_all()

    def phamerated_one(self):
        self.info['phamerated'] = True
        self.info['all'] = False
        self.show_phage_entry()
        self.show_gene_no_entry()
        self.show_starterate_button() 
        self.show_all()    
  
    def unphamerated_one(self):
        self.info['phamerated'] = False
        self.info['all'] = False
        self.show_phage_entry()
        self.show_fasta_entry()
        self.show_gene_no_entry()
        self.show_gene_info_entry()
        self.show_starterate_button()
        self.show_all()
    
    def one_pham(self):
        self.info['phamerated'] = True
        self.info['all'] = False
        self.show_one_pham_entry()
        self.show_starterate_button()
        self.show_all()

    def __init__(self, parent, choice, config):
        Gtk.Dialog.__init__(self, "Starterator", parent.window, 0)
        self.set_border_width(10)
        self.config_info = config
        self.info = {'phage': None,
                    'phamerated' : False,
                    'all' : False,
                    'fasta' : None,
                    'profile': None,
                    'gene_no': None,
                    'start' : None,
                    'stop' : None,
                    'orientation' : None,
                    'pham' : None}
        self.connect("destroy", self.stop_starterator)
        box = self.get_content_area()
        self.box = box
        # self.info_entry_box = Gtk.Box(Gtk.Orientation.VERTICAL, spacing=6)
        # self.box.pack_start(self.info_entry_box, False, False, 0)
        if choice == 'Whole Phamerated Phage':
            self.phamerated_all()
        if choice == 'Whole Unphamerated Phage':
            self.unphamerated_all()
        if choice == 'One Phamerated Gene':
            self.phamerated_one()
        if choice == 'One Unphamerated Gene':
            self.unphamerated_one()
        if choice == 'Pham':
            self.one_pham()


    def db_connect(self):
        try:
            db = database.DB()
            return db
        except:
            dialog = Gtk.MessageDialog(self, 0, Gtk.MessageType.ERROR,
                Gtk.ButtonsType.CANCEL, "Starterator has encountered an error")
            dialog.format_secondary_text("Error connecting to the database. Please check login credentials in Preferences menu")
            dialog.run()
            dialog.destroy()

    def starterate(self, button):
        db = self.db_connect()
        phamgene.check_protein_db(self.config_info['count'])
        phage_name = self.find_phage_in_db(db, str(self.info['phage']))
        print('new phage name', self.info['phage'])
        print('phamerated', self.info['phamerated'])
        print('all', self.info['all'])
        if phage_name == None and self.info["phamerated"] and not self.info["pham"]:
            self.phameratored_exception(self.info["phage"])
            return
        elif self.info["phamerated"]:
            self.info["phage"] = phage_name
        self.progress_label = Gtk.Label('Starterator is starting')
        self.progress_bar = Gtk.ProgressBar()
        self.box.pack_start(self.progress_label, False, False, 0)
        self.box.pack_start(self.progress_bar, False, False, 0)
        self.show_all()
        self.starterate_thread = StarteratorThread(self, db, self.config_info, self.info )
        self.starterate_thread.start()
    
    def stop_starterator(self, button):
        if self.starterate_thread != None:
            self.starterate_thread.cancel()
        self.destroy()
    
    def phameratored_exception(self, name):
        dialog = Gtk.MessageDialog(self, 0, Gtk.MessageType.ERROR,
            Gtk.ButtonsType.OK, "Phage %s could not be found in Phamerator Database!" % name)
        dialog.format_secondary_text("Please check your spelling or choose an Unphameratored option.")
        dialog.run()
        dialog.destroy()
    
    def on_file_clicked(self, button, data):
        name = data[0]
        entry = data[1]
        dialog = Gtk.FileChooserDialog("Please choose the %s file" % name, self,
            Gtk.FileChooserAction.OPEN,
            (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
             Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            print("File selected: " + dialog.get_filename())
            self.info[name] = dialog.get_filename()
            entry.set_text(self.info[name])
            print(name, self.info[name])
        dialog.destroy()

    def on_entry_changed(self, entry, name):
        self.info[name] = entry.get_text()
    
    def on_orientation_toggled(self, button, orientation):
            if button.get_active():
                self.info['orientation'] = orientation
                if orientation == 'F':
                    self.orientation_R.set_active(False)
                elif orientation == 'R':
                    self.orientation_F.set_active(False)
            elif not button.get_active():
                if orientation == 'F':
                    self.orientation_R.set_active(True)
                    self.info['orientation'] = 'R'
                elif orientation == 'R':
                    self.orientation_F.set_active(True)
                    self.info['orientation'] = 'F'


    def get_phage_list(self, db, phage_list):
        cursor = db.cursor()
        cursor.execute('SELECT Name from phage')
        results = cursor.fetchall()
        for row in results:
            phage_list.append(row)

    def find_phage_in_db(self, db, phage):
        results = db.query("SELECT Name\n\
            from phage\n\
            where Name like %s\n\
            or Name like %s \n\
            or Name = %s", (phage+'-%', phage+'_%', phage))
        if len(results) < 1:
            self.info['phamerated'] == False
            return None
        else:
            self.info['phamerated'] == True
            return results[0][0]
    
    def check_pham_number():
        pass

    def update_starterator(self, label_text, progress_amount):
        self.progress_label.set_text(label_text)
        self.progress_bar.set_fraction(progress_amount)

class StarteratorThread(threading.Thread):
    def __init__(self, parent, db, config, info):
        threading.Thread.__init__(self)
        # self.stop = False
        self.parent = parent
        self.config_info = config
        self.info = info
        self.final_file = None
        self.setDaemon(True)
        self.stop_thread = threading.Event()

    def run(self):
        self.starterate()

    def starterate(self):
        try:
            self.final_file = starterate.starterate(self.info, 
                gui=self, event=self.stop_thread)
            # print self.final_file
        except StarteratorError as e:
            self.stop = True
            Gdk.threads_enter()
            dialog = Gtk.MessageDialog(self.parent, 0, Gtk.MessageType.ERROR,
                Gtk.ButtonsType.CANCEL, "Starterator has encountered an error")
            dialog.format_secondary_text(str(e))
            dialog.run()
            dialog.destroy()
            Gdk.threads_leave()
        except:
            if self.stop_thread.is_set():
                return
            print('exception in starterator')
            self.stop = True
            Gdk.threads_enter()
            dialog = Gtk.MessageDialog(self.parent, 0, Gtk.MessageType.ERROR,
                Gtk.ButtonsType.CANCEL, "Startstaerator has encountered an error")
            dialog.format_secondary_text(
            "Please try again.")
            dialog.run()
            dialog.destroy()
            Gdk.threads_leave()
            # self.stop = True
            raise
 
        else:
            if self.stop_thread.is_set():
                # clean up files?
                # show dialog
                return
            print('help after starterator finished')
            Gdk.threads_enter()
            done_dialog = StarteratorFinishedDialog(self.parent, self.final_file)
            response = done_dialog.run()
            if response == Gtk.ResponseType.NO:
                Gtk.main_quit()
            done_dialog.destroy()
            Gdk.threads_leave()
            # self.stop = True

    def cancel(self):
        self.stop_thread.set()

    def update(self, text, amount):
        print('update', text, amount)
        Gdk.threads_enter()
        self.parent.update_starterator(text, amount)
        Gdk.threads_leave()