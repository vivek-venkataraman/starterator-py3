#!/usr/bin/env python

# Copyright (c) 2013, 2014 All Right Reserved, Hatfull Lab, University of Pittsburgh
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
# PARTICULAR PURPOSE.  USE AT YOUR OWN RISK.
#
# Marissa Pacey
# April 4, 2014
# Starting window of Starterator

import MySQLdb
from gi.repository import Gtk
import configparser
import os
import subprocess
import sys
from . import utils
from .uiStarterate import StarteratorEnterInformation
from .uiDialogs import DatabaseInfoDialog, PreferencesDialog
import time
from . import phamgene
from . import phams

"""
GUI Application for Starterator
Process:
    Open Starterator, StarteratorWindow appears
    Choose one of 5 choices, a new dialog pops up based on that choice
    Enter in relevant information, press Starterate
        -Closing this dialog when Starterator is running will cause starterator to stop some time after exit
        -Not instantly, but it checks pretty often-depends on what is happening
    When done, new dialog shows with a link to final report
    If error occurs, dialog shows up
    If unable to connect to the database, DatabaseInfoDialog appears until correct info added
"""

MENU_UI = """
<ui>
  <menubar name='MenuBar'>
    <menu action='File'>
        <menuitem action="ViewReports" />
      <menuitem action='FileQuit' />
    </menu>
    <menu action='Edit'>
      <menuitem action='Preferences' />
    </menu>
    <menu action='Help'>
      <menuitem action='Contents'/>
      <separator />
      <menuitem action='About'/>
    </menu>
  </menubar>
</ui>
"""


class StarteratorWindow:
    def __init__(self):

        self.choices = ['Whole Phamerated Phage', 'Whole Unphamerated Phage', 
                        'One Phamerated Gene', 'One Unphamerated Gene', 'Pham']
  
        # self.set_border_width(10)
        self.config_info = utils.get_config()

        builder = Gtk.Builder()
        builder.add_from_file(utils.glade_file)
        self.window = builder.get_object("StartWindow")
        self.window.set_icon_from_file(utils.icon_file)
        choice_box = builder.get_object('choicebox')
        choice_list = Gtk.ListStore(str)
        for choice in self.choices: choice_list.append([choice])
        choice_box.set_model(choice_list)
        renderer = Gtk.CellRendererText()
        choice_box.pack_start(renderer, True)
        choice_box.add_attribute(renderer, "text", 0)
        choice_box.set_active(0)
        builder.connect_signals(self)
        self.window.show_all()
        # self.vbox.pack_start(vbox, True, True, 0)
        # self.config_info = {}
        # print utils.get_configuration
        self.check_blast_type()
        # db = self.attempt_db_connect()
    
    def on_quit(self, *args, **kwargs):
        Gtk.main_quit()

    def on_view_clicked(self, button):
        print(os.path.abspath(fself.config_info["final_file_dir"]))
        os.system("xdg-open \""+ os.path.abspath(self.config_info["final_file_dir"])+"\"")

    def on_contents_clicked(self, button):
        root_dir = utils.help_files 
        print(root_dir)
        uri = 'ghelp:%s/index.page' % root_dir
        print(uri)
        time_now = int(time.time())

        Gtk.show_uri(None, uri, time_now)

    def on_about_clicked(self, button):
        pass

    def on_pref_clicked(self, data):
        dialog = PreferencesDialog(self, self.config_info)
        response = dialog.run()
        if response == Gtk.ResponseType.APPLY:
            utils.write_to_config_file(self.config_info)
        dialog.destroy()

    def on_fasta_clicked(self, pham_entry):
        pham_no = pham_entry.get_text()
        try:
            pham = phams.Pham(pham_no)
        except:
            dialog = Gtk.MessageDialog(self.window, 0, Gtk.MessageType.ERROR,
                Gtk.ButtonsType.CANCEL, "Pham %s is not a valid pham number" % pham_no)
            dialog.run()
            dialog.destroy()
            return
        dialog = Gtk.FileChooserDialog("Save fasta file for Pham %s" % pham_no, self.window,
            Gtk.FileChooserAction.SAVE,
            (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
             Gtk.STOCK_SAVE, Gtk.ResponseType.OK))
        dialog.set_current_name("Pham%s" %pham_no)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            file_name = dialog.get_filename()
            pham.make_fasta(file_name)
        dialog.destroy()

    def check_blast_type(self):
        blast_dir = self.config_info['blast_dir']
        blast_cmd = self.config_info['blast_dir'] + 'makeblastdb'
        print(os.path.join(blast_dir,'blastp'))
        try:
            print(blast_cmd)
            subprocess.check_call([blast_cmd, '-help'])
            print(blast_cmd, "worked")
            self.config_info['legacy_blast'] = False
        except:
            try:
                print(os.path.join(blast_dir,'formatdb'))
                self.config_info['legacy_blast'] = True
            except:
                dialog = Gtk.MessageDialog(self, 0, Gtk.MessageType.ERROR,
                    Gtk.ButtonsType.CANCEL, "BLAST is not installed")
                dialog.format_secondary_text(
                "Please install BLAST.")
                dialog.run()
                dialog.destroy()
                sys.exit(1)
                
    
    def get_configuration(self):
        config = configparser.RawConfigParser()
        print(utils.STARTERATOR_PATH + "/extras/starterator.config")
        config.read(utils.STARTERATOR_PATH + "/extras/starterator.config")
        print(config)
        self.config_info = dict(config.items('Starterator'))
        # self.config_info['intermediate_file_dir'] = os.path.abspath(self.config_info['intermediate_file_dir'])+ '/'
        # self.config_info['final_file_dir'] = os.path.abspath(self.config_info['final_file_dir']) + '/'
        # self.config_info['protein_db'] = os.path.abspath(self.config_info['protein_db']) + '/'
      

    def show_folders(self, name, readable_name):
        dialog = Gtk.FileChooserDialog("Please choose the folder for storing %s" %readable_name, self,
            Gtk.FileChooserAction.CREATE_FOLDER,
            (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
             Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            print("Open clicked")
            print("File selected: " + dialog.get_filename())
            self.config_info[name] = dialog.get_filename()
            print(name, self.config_info[name])
        elif response == Gtk.ResponseType.CANCEL:
            print("Cancel clicked")
        dialog.destroy()

    def db_connect(self):
        db = MySQLdb.connect(self.config_info['database_server'], 
                self.config_info['database_user'],
                self.config_info['database_password'],
                self.config_info['database_name'])
        return db
    
    def attempt_db_connect(self):
        try:
            print('attempting to connect', self.config_info)
            db = MySQLdb.connect(self.config_info['database_server'], 
                    self.config_info['database_user'],
                    self.config_info['database_password'],
                    self.config_info['database_name'])
            db.close()
        except:
            db_dialog = DatabaseInfoDialog(self, self.config_info['database_server'], 
                    self.config_info['database_name'],
                    self.config_info['database_user'])
            response = db_dialog.run()
            if response == Gtk.ResponseType.OK:
                self.config_info['database_server'] = db_dialog.info['db_server'] 
                self.config_info['database_user'] = db_dialog.info['db_user'] 
                self.config_info['database_password'] = db_dialog.info['db_password']
                self.config_info['database_name'] = db_dialog.info['db_name'] 
            db_dialog.destroy()
            self.attempt_db_connect()
        else:
            utils.write_to_config_file(self.config_info)

    
    def clean_up_files(self):
        for f in os.listdir(self.config_info['intermediate_file_dir']):
            os.remove(os.path.join(self.config_info['intermediate_file_dir'], f))

    def check_files(self):
        if self.config_info["intermediate_file_dir"] == "?":
            self.show_exception("Intermediate files")
            return False
        if self.config_info["final_file_dir"] == "?":
            self.show_exception("Final Reports")
            return False
        # if self.config_info["protein_db"] == "?":
        #     self.show_exception("Protein Database")
        #     return False
        return True
    
    def show_exception(self, name):
        dialog = Gtk.MessageDialog(self, 0, Gtk.MessageType.ERROR,
            Gtk.ButtonsType.OK, "%s folder is not indicated!" % name)
        dialog.format_secondary_text("Please choose folder in the Preferences.")
        dialog.run()
        dialog.destroy()


    def on_choice_changed(self, combo):
        # if self.info_entry_box:
        #     del self.info_entry_box
        if not self.check_files():
            return
        self.config_info['intermediate_file_dir'] = os.path.abspath(self.config_info['intermediate_file_dir'])+ '/'
        self.config_info['final_file_dir'] = os.path.abspath(self.config_info['final_file_dir']) + '/'
        self.config_info['protein_db'] = os.path.abspath(self.config_info['protein_db']) + '/'
        tree_iter = combo.get_active_iter()
        if tree_iter != None:
            model = combo.get_model()
            choice = model[tree_iter][0]
            dialog = StarteratorEnterInformation(self, choice, self.config_info)