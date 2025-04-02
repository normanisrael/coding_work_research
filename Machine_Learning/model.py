import numpy as np
import os, shutil
from keras import layers
from keras import models
from keras import optimizers
from keras.preprocessing.image import ImageDataGenerator

TRAIN_DIR = '/codes/data/train'
VAL_DIR = '/codes/data/validation'

input_shape = (150,150,3)
model = models.Sequential()
model.add(layers.Conv2D(32, (3,3), activation='relu', input_shape=input_shape))
model.add(layers.MaxPooling2D((2,2)))
model.add(layers.Dropout(0.25))
model.add(layers.Conv2D(64, (3,3), activation='relu'))
model.add(layers.MaxPooling2D((2,2)))
model.add(layers.Dropout(0.25))
model.add(layers.Conv2D(128, (3,3), activation='relu'))
model.add(layers.MaxPooling2D((2,2)))
model.add(layers.Dropout(0.25))
model.add(layers.Conv2D(128, (3,3), activation='relu'))
model.add(layers.MaxPooling2D((2,2)))
model.add(layers.Dropout(0.25))
model.add(layers.Flatten())
model.add(layers.Dense(512, activation='relu'))
model.add(layers.Dropout(0.5)) #0.8 for diagnostic 0
model.add(layers.Dense(1, activation='sigmoid')) #sigmoid initially

model.compile(loss='binary_crossentropy', optimizer=optimizers.RMSprop(lr=1e-3), metrics=['acc'])

train_gen_rescale = ImageDataGenerator(rescale=1./255,
                                        width_shift_range=0.2,
                                        height_shift_range=0.2,
                                        shear_range=0.2,
                                        zoom_range=0.2,
                                        horizontal_flip=True,
                                        fill_mode='nearest')

val_gen_rescale = ImageDataGenerator(rescale=1./255)

train_gen = train_gen_rescale.flow_from_directory(TRAIN_DIR,
                                              target_size=(150,150),
                                              batch_size=15,
                                              class_mode='binary')

val_gen = val_gen_rescale.flow_from_directory(VAL_DIR,
                                             target_size=(150,150),
                                             batch_size=15,
                                             class_mode='binary')

history = model.fit_generator(train_gen,
                             steps_per_epoch=100,
                             epochs=100,
                             validation_data=val_gen,
                             validation_steps=3)

training_accuracy = history.history['acc']
validation_accuracy = history.history['val_acc']
training_loss = history.history['loss']
validation_loss = history.history['val_loss']

epochs = range(1, len(training_accuracy) + 1)

model_diagnostics = open("model_diagnostics_0.dat", "a")
#model_diagnostics.write('{:<10s}{:10s}{:10s}{:15s}{:>15s} \n'.format('epoch', 'training accuracy', 'validation accuracy', 'training loss', 'validation loss'))
for i in range(0, len(epochs)):
    model_diagnostics.write('{:<3d}{:10f}{:10f}{:15.9f}{:>15.9f} \n'.format(epochs[i], training_accuracy[i], validation_accuracy[i], training_loss[i], validation_loss[i]))

model_diagnostics.close()

  #print('{:<3d}{:10f}{:10f}{:15.9f}{:>15.9f}'.format(epochs[i], acc[i], val_acc[i], loss[i], val_loss[i]))
